# 二倍体短读段定相服务（Diploid Short-Read Phasing）

基于 **Python + FastAPI + OR-Tools (CP-SAT)** 的二倍体短读段定相服务。
将每条读段归入两条互补单倍型之一，以**不符观测的权重之和**为代价，
用整数规划求**已证明的全局最小值**（MEC，Minimum Error Correction）。

- 按共覆盖关系把位点分成若干连通块，逐块独立求解，**不推断块间相位**；
- 每块首位的单倍型等位固定为 `0`，消除互补对称；
- 返回每块的最小代价、按位点顺序**字典序最小**的最优单倍型，
  以及在所有最优解中**可翻转**（两种等位都出现在某个最优解中）的位点；
- 求解器未能在时限内证明最优时返回 `undecided`，**绝不给出伪最优结果**；
- 重复位点、非法读段等请求错误返回 **422**。

## 快速启动（Docker Compose）

```bash
docker compose up --build
```

服务监听 <http://localhost:8000>：

- `POST /phase` — 定相接口
- `GET /health` — 健康检查
- <http://localhost:8000/docs> — Swagger UI 交互式文档（可直接试调）

## 本地运行（不用 Docker）

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/uvicorn app.main:app --host 0.0.0.0 --port 8000
```

运行测试：

```bash
.venv/bin/pip install -r requirements-dev.txt
.venv/bin/pytest tests/ -q
```

## API 说明

### `POST /phase`

**请求体**

```json
{
  "sites": [
    {"id": "s1", "position": 100}
  ],
  "reads": [
    {
      "id": "r1",
      "observations": [
        {"site": "s1", "allele": 0, "weight": 30}
      ]
    }
  ]
}
```

| 字段 | 约束 |
| --- | --- |
| `sites` | 至少 1 个；`id` 唯一；`position` 严格递增（即按位置排序且无重复） |
| `reads` | 可为空；每条读段至少 2 条观测 |
| `observations[].site` | 必须出现在 `sites` 中；同一读段内不得重复 |
| `observations[].allele` | 整数 `0` 或 `1` |
| `observations[].weight` | 整数 `1`～`60` |
| `reads[].id` | 可选，仅供客户端追踪 |

违反以上约束一律返回 **422**。

**响应（已证明最优）** — HTTP 200

| 字段 | 含义 |
| --- | --- |
| `status` | `"optimal"` |
| `total_cost` | 各块最小代价之和 |
| `blocks[].site_ids` | 块内位点，顺序与请求一致 |
| `blocks[].min_cost` | 该块最小代价（不符观测的权重之和） |
| `blocks[].haplotype` | 字典序最小的最优单倍型（0/1 字符串，与 `site_ids` 对齐，首位固定为 `0`） |
| `blocks[].flippable_site_ids` | 在所有最优解中可取另一等位的位点（块首位因对称性破缺固定，永不包含） |

**响应（未决）** — HTTP 200，`blocks` 为 `null`，不给出任何伪最优结果：

```json
{
  "status": "undecided",
  "total_cost": null,
  "blocks": null,
  "detail": "primary solve did not prove optimality (status=FEASIBLE)"
}
```

## 实际请求与响应示例

```bash
curl -s -X POST http://localhost:8000/phase \
  -H 'Content-Type: application/json' \
  -d '{
    "sites": [
      {"id": "s1", "position": 100},
      {"id": "s2", "position": 200},
      {"id": "s3", "position": 300},
      {"id": "s4", "position": 400},
      {"id": "s5", "position": 500}
    ],
    "reads": [
      {"id": "r1", "observations": [{"site": "s1", "allele": 0, "weight": 30}, {"site": "s2", "allele": 0, "weight": 30}]},
      {"id": "r2", "observations": [{"site": "s1", "allele": 0, "weight": 10}, {"site": "s2", "allele": 1, "weight": 10}]},
      {"id": "r3", "observations": [{"site": "s2", "allele": 0, "weight": 30}, {"site": "s3", "allele": 1, "weight": 30}]},
      {"id": "r4", "observations": [{"site": "s4", "allele": 1, "weight": 20}, {"site": "s5", "allele": 1, "weight": 20}]},
      {"id": "r5", "observations": [{"site": "s4", "allele": 1, "weight": 20}, {"site": "s5", "allele": 0, "weight": 20}]}
    ]
  }'
```

实际响应（`s1–s3` 与 `s4–s5` 互不共覆盖，分为两块；块 2 中 `r4`/`r5`
等权冲突，`s5` 在另一最优解中可取 `1`，故可翻转）：

```json
{
  "status": "optimal",
  "total_cost": 30,
  "blocks": [
    {
      "index": 0,
      "site_ids": ["s1", "s2", "s3"],
      "min_cost": 10,
      "haplotype": "001",
      "flippable_site_ids": []
    },
    {
      "index": 1,
      "site_ids": ["s4", "s5"],
      "min_cost": 20,
      "haplotype": "00",
      "flippable_site_ids": ["s5"]
    }
  ],
  "detail": null
}
```

非法请求示例（位点 `id` 重复）：

```bash
curl -s -X POST http://localhost:8000/phase \
  -H 'Content-Type: application/json' \
  -d '{"sites": [{"id": "s1", "position": 100}, {"id": "s1", "position": 200}], "reads": []}'
```

返回 HTTP **422**：

```json
{
  "detail": [
    {
      "type": "value_error",
      "loc": ["body"],
      "msg": "Value error, duplicate site id(s): ['s1']"
    }
  ]
}
```

## 求解与配置

数学模型：对每个位点 `i` 设 0/1 变量 `h_i`（单倍型 0 的等位，单倍型 1 为其互补），
对每条读段 `r` 设 0/1 变量 `x_r`（归入哪条单倍型）；观测 `(r, i)` 的等位
`a ≠ h_i XOR x_r` 时产生代价 `weight`，最小化总代价。每块固定 `h_0 = 0`
消除互补对称。字典序最小单倍型与可翻转位点通过“代价优先 + 次级目标”的
组合目标函数在已证明的最优解集合内精确求得。

环境变量（可在 `docker-compose.yml` 中调整）：

| 变量 | 默认值 | 含义 |
| --- | --- | --- |
| `PHASE_MAX_SOLVER_SECONDS` | `10` | 单次 CP-SAT 求解的时间上限（秒） |
| `PHASE_MAX_TOTAL_SECONDS` | `120` | 单个请求的总时间预算（秒） |
| `PHASE_SOLVER_WORKERS` | CPU 核数 | CP-SAT 搜索线程数 |

MEC 为 NP-hard 问题。参考规模（12 核、默认配置）：约 150 个位点、1200 条
读段的单块实例可在总预算内完成全部证明；更大或更嘈杂的实例可能返回
`undecided`，此时可调高上述时限。

## 项目结构

```
app/
  main.py      # FastAPI 入口，POST /phase
  schemas.py   # 请求/响应模型与 422 校验
  phaser.py    # 连通块分解 + CP-SAT 精确求解
tests/
  test_phase.py  # 422 校验、正确性（含暴力枚举对照）、undecided 行为
Dockerfile
docker-compose.yml
requirements.txt / requirements-dev.txt
```
