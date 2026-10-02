"""核心定相算法：MEC（最小错误校正）整数规划 + 连通块分解。

对每条读段 r 引入 0/1 变量 x_r（归入哪条单倍型），对每个位点 i 引入
0/1 变量 h_i（单倍型 0 在该位点的等位，单倍型 1 为其互补）。观测 (r, i)
的等位 a 与 h_i XOR x_r 不一致时产生代价 weight。用 OR-Tools CP-SAT
求全局最小值；只有状态为 OPTIMAL（已证明最优）才返回结果，否则抛
UndecidedError，由上层返回 "undecided"。

次要分析全部采用“代价优先 + 次级目标”的组合目标函数
（cost * M + secondary，M 严格大于次级目标上界），保证每次求解都在
全局最优解集合内做精确优化，且界传播与主求解一样强：

- 字典序最小单倍型：按 32 位分块，最小化 cost * 2^32 + 块内二进制数，
  已确定的块用硬约束固定；
- 可翻转位点：以主求解得到的最优解为参照，每轮最小化
  cost * (n+1) -（未确认位点的翻转数），即在不牺牲最优性的前提下
  尽可能多地翻转未确认位点；一轮翻转数为 0 时，剩余位点均不可翻转。
"""

from __future__ import annotations

import os
import time
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

from ortools.sat.python import cp_model

from .schemas import BlockResult, PhaseRequest

# 字典序最小化时每个二进制块的位数（对应组合目标中代价的乘数 2^32）
_LEX_CHUNK_BITS = 32


class UndecidedError(Exception):
    """求解器未能在时限内证明最优（或完成等价的最优性分析）。"""


@dataclass(frozen=True)
class SolverConfig:
    max_solver_seconds: float = 10.0   # 单次 CP-SAT 求解的时间上限
    max_total_seconds: float = 120.0   # 单个请求的总时间预算
    num_workers: int = 0               # CP-SAT 搜索线程数，0 表示取 CPU 核数


def config_from_env() -> SolverConfig:
    def _f(name: str, default: float) -> float:
        try:
            return float(os.environ.get(name, default))
        except ValueError:
            return default

    return SolverConfig(
        max_solver_seconds=_f("PHASE_MAX_SOLVER_SECONDS", 10.0),
        max_total_seconds=_f("PHASE_MAX_TOTAL_SECONDS", 120.0),
        num_workers=int(_f("PHASE_SOLVER_WORKERS", 0)),
    )


# ---------------------------------------------------------------------------
# 连通块分解
# ---------------------------------------------------------------------------

def _connected_components(
    n_sites: int, reads_sites: Sequence[Sequence[int]]
) -> List[List[int]]:
    """按共覆盖关系对位点做并查集连通块分解，返回按块内最小位点排序的块列表。"""
    parent = list(range(n_sites))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for sites in reads_sites:
        for s in sites[1:]:
            union(sites[0], s)

    blocks: Dict[int, List[int]] = {}
    for i in range(n_sites):
        blocks.setdefault(find(i), []).append(i)
    return sorted(blocks.values(), key=lambda b: b[0])


# ---------------------------------------------------------------------------
# 单块求解
# ---------------------------------------------------------------------------

# 读段观测展开后的形式：(局部位点下标, 等位, 权重)
_Obs = Tuple[int, int, int]


def _phase_block(
    n: int,
    block_reads: Sequence[Sequence[_Obs]],
    config: SolverConfig,
    deadline: float,
) -> Tuple[int, List[int], List[int]]:
    """求解一个连通块。

    返回 (最小代价, 字典序最小的最优单倍型, 可翻转位点的局部下标列表)。
    单倍型首位固定为 0 以消除互补对称。
    """
    if n == 1:
        # 孤立位点：无任何观测，代价为 0，首位固定为 0 且不可翻转
        return 0, [0], []

    model = cp_model.CpModel()
    h = [model.new_bool_var(f"h_{i}") for i in range(n)]
    model.add(h[0] == 0)  # 对称性破缺：块内首位固定为 0

    cost_terms = []
    for r_idx, obs_list in enumerate(block_reads):
        x = model.new_bool_var(f"x_{r_idx}")  # 0: 归入单倍型0, 1: 归入单倍型1
        for site, allele, weight in obs_list:
            e = model.new_bool_var(f"e_{r_idx}_{site}")
            # 读段 r 在 site 上被解释的等位为 h[site] XOR x。
            # 只需下界约束：目标最小化会把 e 压到最小。
            if allele == 0:  # e >= h XOR x
                model.add(e >= h[site] - x)
                model.add(e >= x - h[site])
            else:            # e >= 1 - (h XOR x)，即 h == x 时 e = 1
                model.add(e >= h[site] + x - 1)
                model.add(e >= 1 - h[site] - x)
            cost_terms.append(weight * e)

    total_cost = sum(cost_terms) if cost_terms else 0
    model.minimize(total_cost)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = config.max_solver_seconds
    solver.parameters.num_search_workers = config.num_workers or (os.cpu_count() or 1)

    def _solve(what: str) -> None:
        """求解当前模型，要求必须证明最优，否则抛 UndecidedError。"""
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            raise UndecidedError("request time budget exceeded")
        solver.parameters.max_time_in_seconds = min(
            config.max_solver_seconds, remaining
        )
        st = solver.solve(model)
        if st != cp_model.OPTIMAL:
            raise UndecidedError(
                f"{what} did not prove optimality (status={solver.status_name(st)})"
            )

    _solve("primary solve")
    optimal_cost = int(round(solver.objective_value))
    # 参照最优解：可翻转性定义为“存在取另一等位的全局最优解”，与参照无关
    reference = [int(round(solver.value(h[i]))) for i in range(n)]

    # --- 可翻转位点 -----------------------------------------------------------
    # 每轮在未确认的位点上最大化翻转数（代价优先），翻转到的位点即被实际
    # 最优解见证为可翻转；一轮无任何翻转时，剩余位点在所有最优解中恒定。
    flippable = set()
    undecided_sites = set(range(1, n))
    while undecided_sites:
        flip_count = sum(
            h[i] if reference[i] == 0 else 1 - h[i] for i in undecided_sites
        )
        model.clear_objective()
        model.minimize(total_cost * (n + 1) - flip_count)
        _solve("flippability solve")
        gained = {
            i
            for i in undecided_sites
            if int(round(solver.value(h[i]))) != reference[i]
        }
        if not gained:
            break
        flippable |= gained
        undecided_sites -= gained
    model.clear_objective()

    # --- 字典序最小的最优单倍型 ---------------------------------------------
    # 按 32 位分块，最小化 cost * 2^32 + 块内位点构成的二进制数；已确定的
    # 块用硬约束固定（硬约束不阻碍 presolve，效果远好于 assumptions）。
    haplotype = [0] * n
    pos = 1
    while pos < n:
        end = min(n, pos + _LEX_CHUNK_BITS)
        chunk_value = sum(h[i] * (1 << (end - 1 - i)) for i in range(pos, end))
        model.clear_objective()
        model.minimize(total_cost * (1 << _LEX_CHUNK_BITS) + chunk_value)
        _solve("lexicographic solve")
        value = int(round(solver.objective_value)) - optimal_cost * (1 << _LEX_CHUNK_BITS)
        for i in range(pos, end):
            haplotype[i] = (value >> (end - 1 - i)) & 1
        model.add(chunk_value == value)
        pos = end

    return optimal_cost, haplotype, sorted(flippable)


# ---------------------------------------------------------------------------
# 入口
# ---------------------------------------------------------------------------

def run_phasing(request: PhaseRequest, config: SolverConfig) -> List[BlockResult]:
    """对请求做连通块分解并逐块求解，返回各块结果。"""
    site_index = {s.id: i for i, s in enumerate(request.sites)}
    reads_sites: List[List[int]] = [
        [site_index[o.site] for o in read.observations] for read in request.reads
    ]

    blocks = _connected_components(len(request.sites), reads_sites)

    # 把读段按块归组（一条读段覆盖的所有位点必在同一块内）
    block_of_site: Dict[int, int] = {}
    for b_idx, block in enumerate(blocks):
        for s in block:
            block_of_site[s] = b_idx
    block_reads: List[List[_Obs]] = [[] for _ in blocks]
    for read, sites in zip(request.reads, reads_sites):
        b_idx = block_of_site[sites[0]]
        local = {g: l for l, g in enumerate(blocks[b_idx])}
        block_reads[b_idx].append(
            [(local[g], obs.allele, obs.weight) for g, obs in zip(sites, read.observations)]
        )

    deadline = time.monotonic() + config.max_total_seconds
    results: List[BlockResult] = []
    for b_idx, block in enumerate(blocks):
        cost, haplotype, flippable = _phase_block(
            len(block), block_reads[b_idx], config, deadline
        )
        site_ids = [request.sites[g].id for g in block]
        results.append(
            BlockResult(
                index=b_idx,
                site_ids=site_ids,
                min_cost=cost,
                haplotype="".join(map(str, haplotype)),
                flippable_site_ids=[site_ids[i] for i in flippable],
            )
        )
    return results
