# 二倍体短读段定相服务（Diploid Short-Read Phasing）

基于 **Python + FastAPI + OR-Tools CP-SAT** 的二倍体短读段（short reads）定相服务，
提供一个可直接调用的 `POST /phase` HTTP API。使用 Docker Compose 一键启动。

## 解决什么问题

给定一组**按位置排序、唯一的杂合双等位位点**以及一组短读段：

- 每条读段覆盖**至少两个已知位点**，读段在每个覆盖位点上有一个 `0/1` 等位观测
  和一个 **1～60 的整数权重**（可表示测序质量/置信度）；
- 每条读段必须被**整体**归入两条互补单倍型之一；
- 观测与归属不符时，代价加上该观测的权重；
- 目标是**全局最小化不符观测的权重之和**（加权 MEC, Minimum Error Correction）。

本服务用 OR-Tools CP-SAT 对该问题做**精确求解**（枚举/证明最优，而非启发式近似）。

### 连通块与互补对称

- 按读段的**共覆盖关系**（两条位点被同一条读段覆盖即在同一块）把位点划分成连通块；
- 不同块之间没有任何读段联系，**不推断块间相位**，各块独立求解、独立返回；
- 把每个块按位置排序后的**首位固定为 0**，消除“所有位点 + 所有读段归属整体取反”
  的互补对称，保证解的表示唯一。

### 返回内容

每个连通块返回：

1. **`min_cost`**：该块被证明的最小代价；
2. **`haplotype`**：在所有最优解中，按位点顺序**字典序最小**的最优单倍型（0/1 序列）；
3. **`flippable_sites`**：在**所有最优解**中可以翻转的位点。某位点同时存在最优解
   使其取 0 与取 1 时即为可翻转（不要求其余位点保持不变；必须成对翻转的位点也
   会被列出）。锚点首位恒为 0，永不可翻转。

### 未决（undecided）语义

求解器在给定时间内**未能证明最优性**时，服务**绝不返回伪最优结果**：

- 对应连通块 `"status": "undecided"`，且 `min_cost` / `haplotype` /
  `flippable_sites` 全部为 `null`，只返回块内位点下标；
- 顶层 `status` 只要存在任何未决块即为 `"undecided"`；
- 已证明最优的块仍正常给出完整结果。

## 一键启动

需要本机已安装 Docker 与 Docker Compose（`docker compose` 子命令）。

```bash
docker compose up --build -d
```

服务监听 **http://localhost:8000** 。停止：

```bash
docker compose down
```

健康检查：

```bash
curl http://localhost:8000/health
# {"status":"ok"}
```

交互式 API 文档（Swagger UI）：浏览器打开 <http://localhost:8000/docs>。

## API

### `POST /phase`

**请求体字段**

| 字段 | 说明 |
| --- | --- |
| `sites` | 非空位点数组，元素为 `{"position": <非负整数>}`；必须按 `position` **严格递增且无重复** |
| `reads` | 非空读段数组 |
| `reads[].id` | 读段唯一字符串 ID（1～128 字符），不可重复 |
| `reads[].observations` | 至少 **2** 个观测；同一条读段内 `site` 不可重复 |
| `observations[].site` | 位点在 `sites` 数组中的下标（从 0 开始），不可越界 |
| `observations[].allele` | 观测等位，只能是 `0` 或 `1` |
| `observations[].weight` | 整数权重，**1～60** |
| `time_limit_seconds` | 可选，每次求解器调用的最长秒数（1～120，默认 20） |

**错误处理**：位点未排序/重复、读段 ID 重复、读段覆盖位点数不足、观测位点越界、
同读段重复位点、等位非 0/1、权重超出 1～60、多出未知字段等，一律返回 **HTTP 422**，
响应体为 FastAPI 标准的 `detail` 校验错误列表。

### 示例 1：两个连通块，块间不定相

5 个位点。`r1/r2` 把位点 0、1、2 连成一块（观测自洽，零代价，单倍型 `010`），
`r3` 单独覆盖位点 3、4 形成第二块（单倍型 `01`）。两块之间没有读段联系，
因此**它们的相对相位不做任何推断**。

请求：

```bash
curl -s -X POST http://localhost:8000/phase \
  -H 'Content-Type: application/json' \
  -d '{
    "sites": [{"position": 100}, {"position": 200}, {"position": 300},
              {"position": 400}, {"position": 500}],
    "reads": [
      {"id": "r1", "observations": [
        {"site": 0, "allele": 0, "weight": 30},
        {"site": 1, "allele": 1, "weight": 30}]},
      {"id": "r2", "observations": [
        {"site": 1, "allele": 1, "weight": 20},
        {"site": 2, "allele": 0, "weight": 20}]},
      {"id": "r3", "observations": [
        {"site": 3, "allele": 0, "weight": 10},
        {"site": 4, "allele": 1, "weight": 10}]}
    ]
  }'
```

实际响应：

```json
{
    "status": "optimal",
    "components": [
        {
            "status": "optimal",
            "sites": [0, 1, 2],
            "min_cost": 0,
            "haplotype": [0, 1, 0],
            "flippable_sites": []
        },
        {
            "status": "optimal",
            "sites": [3, 4],
            "min_cost": 0,
            "haplotype": [0, 1],
            "flippable_sites": []
        }
    ]
}
```

### 示例 2：冲突代价、字典序最小解与可翻转位点

3 个位点上有 5 条权重均为 1 的读段。`r1/r4` 在位点 1 上矛盾，`r2/r5` 在位点 2
上矛盾：选择单倍型 `[0,0,0]` 或 `[0,1,1]` 都要违背 2 个权重，最小代价均为 2。
因此字典序最小解是 `[0,0,0]`，位点 1、2 在最优解之间成对翻转，均列入
`flippable_sites`（位点 0 是锚点，不可翻转）。

请求：

```bash
curl -s -X POST http://localhost:8000/phase \
  -H 'Content-Type: application/json' \
  -d '{
    "sites": [{"position": 10}, {"position": 20}, {"position": 30}],
    "reads": [
      {"id": "r1", "observations": [{"site": 0, "allele": 0, "weight": 1},
                                    {"site": 1, "allele": 1, "weight": 1}]},
      {"id": "r2", "observations": [{"site": 0, "allele": 0, "weight": 1},
                                    {"site": 2, "allele": 1, "weight": 1}]},
      {"id": "r3", "observations": [{"site": 1, "allele": 0, "weight": 1},
                                    {"site": 2, "allele": 0, "weight": 1}]},
      {"id": "r4", "observations": [{"site": 0, "allele": 0, "weight": 1},
                                    {"site": 1, "allele": 0, "weight": 1}]},
      {"id": "r5", "observations": [{"site": 0, "allele": 0, "weight": 1},
                                    {"site": 2, "allele": 0, "weight": 1}]}
    ]
  }'
```

实际响应：

```json
{
    "status": "optimal",
    "components": [
        {
            "status": "optimal",
            "sites": [0, 1, 2],
            "min_cost": 2,
            "haplotype": [0, 0, 0],
            "flippable_sites": [1, 2]
        }
    ]
}
```

### 示例 3：权重裁决

权重 60 的深读段与权重 1 的浅读段矛盾时，最优解违背浅读段，代价 1：

```json
// 请求
{
  "sites": [{"position": 100}, {"position": 200}],
  "reads": [
    {"id": "deep", "observations": [{"site": 0, "allele": 0, "weight": 60},
                                    {"site": 1, "allele": 1, "weight": 60}]},
    {"id": "shallow", "observations": [{"site": 0, "allele": 0, "weight": 1},
                                       {"site": 1, "allele": 0, "weight": 1}]}
  ]
}
```

```json
// 实际响应
{
    "status": "optimal",
    "components": [
        {
            "status": "optimal",
            "sites": [0, 1],
            "min_cost": 1,
            "haplotype": [0, 1],
            "flippable_sites": []
        }
    ]
}
```

### 示例 4：非法输入返回 422

位点位置重复（也适用于未排序、读段覆盖位点不足、等位/权重越界、位点越界等）：

```bash
curl -i -X POST http://localhost:8000/phase \
  -H 'Content-Type: application/json' \
  -d '{
    "sites": [{"position": 100}, {"position": 100}],
    "reads": [{"id": "r1", "observations": [
      {"site": 0, "allele": 0, "weight": 1},
      {"site": 1, "allele": 1, "weight": 1}]}]
  }'
```

```
HTTP/1.1 422 Unprocessable Entity
```

```json
{
  "detail": [
    {
      "type": "value_error",
      "loc": ["body"],
      "msg": "Value error, sites 必须按 position 严格递增且不得包含重复位置",
      "input": { "...": "..." }
    }
  ]
}
```

## 本地开发（不使用 Docker）

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload --port 8000
```

运行测试（包含与暴力枚举的随机对拍）：

```bash
pip install pytest httpx
pytest -q
```

## 项目结构

```
.
├── app/
│   ├── __init__.py
│   ├── main.py       # FastAPI 路由 / 请求编排
│   ├── models.py     # Pydantic 请求/响应模型与输入校验（422）
│   └── solver.py     # 连通块划分 + CP-SAT 精确求解（代价/字典序/可翻转）
├── tests/
│   └── test_solver.py
├── Dockerfile
├── docker-compose.yml
└── requirements.txt
```

## 求解方法摘要

对每个连通块建立 CP-SAT 布尔模型：

- `x_i ∈ {0,1}`：位点 `i` 在单倍型 0 上的等位；
- `y_r ∈ {0,1}`：读段 `r` 归属的单倍型；
- 不符指示 `m = x_i XOR allele XOR y_r`，目标
  `min Σ m · weight`；
- 约束 `x_0 = 0` 固定块内首位。

每个块依次执行：① 最小化代价（要求 **OPTIMAL** 状态，否则该块未决）；
② 从左到右逐位在“代价 = 最优值”的可行解中尝试取 0，求字典序最小最优单倍型；
③ 对每个非锚点位点测试“强制取字典序解的反值后是否仍存在最优解”，得到
可翻转位点。任一次调用超时/状态未知都会使对应块返回 `undecided`，
不会输出未经证明的数值。
