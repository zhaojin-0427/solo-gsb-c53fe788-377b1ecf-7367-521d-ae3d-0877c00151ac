"""二倍体短读段定相核心逻辑。

问题形式化（加权 MEC, Minimum Error Correction）：
    - 每条读段 r 被整体归入两条互补单倍型之一 y_r ∈ {0, 1}；
    - 每个位点 i 的单倍型 0 等位取值为 x_i ∈ {0, 1}；
    - 读段 r 在位点 i 上观测到等位 a（权重 w）与归属相符当且仅当
          一致条件: x_i XOR (1 - a) XOR y_r == 0
      即不符指示 m = x_i XOR a XOR y_r；
    - 目标是最小化 sum(m * w)（不符观测权重之和）。

按读段的共覆盖关系把位点分成连通块，块之间没有任何共同读段，
因此块间相位不可推断，各自独立求解。每块按位置排序后的首位固定为 0，
消除 (x, y) 整体取反带来的互补对称。

可翻转位点定义：在保持锚点首位为 0 的全部最优解中，该位点既可以取 0
也可以取 1（不要求其余位点不变；成对翻转的位点也算可翻转）。
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from ortools.sat.python import cp_model

from .models import Read

# 默认每次求解器调用的墙上时间上限（秒）。
DEFAULT_TIME_LIMIT = 20.0


@dataclass
class SolvedComponent:
    sites: List[int]                        # 块内位点的全局下标（按位置升序）
    min_cost: Optional[int] = None          # 未决时为 None
    haplotype: Optional[List[int]] = None   # 字典序最小的最优单倍型
    flippable_sites: Optional[List[int]] = None
    decided: bool = True


# --------------------------------------------------------------------------- #
# 连通块划分
# --------------------------------------------------------------------------- #
class DSU:
    def __init__(self, n: int) -> None:
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


def split_components(
    n_sites: int, reads: List[Read]
) -> List[Tuple[List[int], List[Read]]]:
    """按读段共覆盖关系划分连通块。

    每个元素为 (块内位点全局下标升序, 覆盖该块的读段)；
    块之间按块内最小位点升序排列。没有读段覆盖的孤立位点各自成块。
    """
    dsu = DSU(n_sites)
    for read in reads:
        covered = [obs.site for obs in read.observations]
        first = covered[0]
        for site in covered[1:]:
            dsu.union(first, site)

    groups: Dict[int, List[int]] = {}
    for site in range(n_sites):
        groups.setdefault(dsu.find(site), []).append(site)

    ordered_roots = sorted(groups, key=lambda root: min(groups[root]))

    result: List[Tuple[List[int], List[Read]]] = []
    for root in ordered_roots:
        comp_sites = sorted(groups[root])
        site_set = set(comp_sites)
        # 读段覆盖的所有位点一定在同一个块内（union 保证），取首个观测定位即可。
        comp_reads = [r for r in reads if r.observations[0].site in site_set]
        result.append((comp_sites, comp_reads))
    return result


# --------------------------------------------------------------------------- #
# CP-SAT 建模
# --------------------------------------------------------------------------- #
@dataclass
class ComponentModel:
    model: cp_model.CpModel
    x: List[cp_model.IntVar]   # 位点变量（块内下标）
    y: List[cp_model.IntVar]   # 读段归属变量
    cost: cp_model.IntVar


def build_component_model(
    comp_sites: List[int],
    comp_reads: List[Read],
    fixed: Optional[Dict[int, int]] = None,
    force_cost: Optional[int] = None,
    minimize: bool = False,
) -> ComponentModel:
    """构建一个连通块的 CP-SAT 模型。

    fixed:      额外固定的位点变量（块内下标 -> 0/1），锚点首位恒固定为 0；
    force_cost: 要求总代价恰好等于该值（用于最优解集合上的可行性测试）；
    minimize:   是否最小化总代价（仅首轮求最优代价时使用）。
    """
    fixed = fixed or {}
    n = len(comp_sites)
    local_index = {global_idx: local for local, global_idx in enumerate(comp_sites)}

    model = cp_model.CpModel()

    # 恒为 1 的常量，供 AddBoolXOr 表达等式型 XOR 约束。
    one = model.NewBoolVar("one")
    model.Add(one == 1)

    x = [model.NewBoolVar(f"x_{comp_sites[i]}") for i in range(n)]
    y = [model.NewBoolVar(f"y_{r.id}") for r in comp_reads]

    # 消除互补对称：按位置排序后的首位固定为 0。
    model.Add(x[0] == 0)
    for local, value in fixed.items():
        model.Add(x[local] == value)

    mismatch_terms: List[cp_model.LinearExpr] = []
    for r_idx, read in enumerate(comp_reads):
        for obs in read.observations:
            li = local_index[obs.site]
            # 不符指示 m 与 x、等位、读段归属满足 m = x XOR allele XOR y。
            m = model.NewBoolVar(f"m_{read.id}_{obs.site}")
            # x XOR allele 的文字：allele=0 时为 x，allele=1 时为 not x。
            allele_lit = x[li] if obs.allele == 0 else x[li].Not()
            # AddBoolXOr 要求奇数个文字为真；加入常量文字 1 后，
            # m XOR (x XOR allele) XOR y = 1 等价于 m = x XOR allele XOR y。
            model.AddBoolXOr([m, allele_lit, y[r_idx], one])
            mismatch_terms.append(obs.weight * m)

    total_weight = sum(obs.weight for read in comp_reads for obs in read.observations)
    cost = model.NewIntVar(0, total_weight, "cost")
    if mismatch_terms:
        model.Add(cost == sum(mismatch_terms))
    else:
        model.Add(cost == 0)

    if force_cost is not None:
        model.Add(cost == force_cost)
    if minimize:
        model.Minimize(cost)

    return ComponentModel(model=model, x=x, y=y, cost=cost)


def run_solver(
    cm: ComponentModel,
    time_limit: float,
    hint_x: Optional[List[int]] = None,
) -> Tuple[str, Optional[cp_model.CpSolver]]:
    """求解并返回 (状态类别, solver)。

    状态类别：
      "feasible"   —— OPTIMAL 或 FEASIBLE（可行赋值可用）；
      "infeasible" —— INFEASIBLE；
      "unknown"    —— UNKNOWN / MODEL_INVALID 等无法下结论的情形。
    """
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = 8
    if hint_x is not None:
        for var, value in zip(cm.x, hint_x):
            cm.model.AddHint(var, value)

    status = solver.Solve(cm.model)
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return "feasible", solver
    if status == cp_model.INFEASIBLE:
        return "infeasible", solver
    return "unknown", None


# --------------------------------------------------------------------------- #
# 单块求解
# --------------------------------------------------------------------------- #
def solve_component(
    comp_sites: List[int],
    comp_reads: List[Read],
    time_limit: float,
) -> SolvedComponent:
    n = len(comp_sites)

    # 孤立位点（没有任何读段覆盖）：锚点定为 0，代价 0，无翻转自由度。
    if not comp_reads:
        return SolvedComponent(
            sites=comp_sites,
            min_cost=0,
            haplotype=[0] * n,
            flippable_sites=[],
            decided=True,
        )

    # 第 1 步：求最小代价，必须由求解器证明最优。
    cm = build_component_model(comp_sites, comp_reads, minimize=True)
    kind, solver = run_solver(cm, time_limit)
    if kind != "feasible" or solver.StatusName() != "OPTIMAL":
        return SolvedComponent(sites=comp_sites, decided=False)
    opt_cost = solver.Value(cm.cost)
    x_values = [solver.Value(v) for v in cm.x]

    # 第 2 步：在所有最优解中求字典序最小的单倍型。
    # 从左到右逐位尝试固定为 0；不可行则该位必为 1；超时即未决。
    fixed: Dict[int, int] = {}
    for i in range(1, n):
        trial = dict(fixed)
        trial[i] = 0
        cm = build_component_model(comp_sites, comp_reads, fixed=trial, force_cost=opt_cost)
        kind, solver = run_solver(cm, time_limit, hint_x=x_values)
        if kind == "feasible":
            fixed[i] = 0
            x_values = [solver.Value(v) for v in cm.x]
        elif kind == "infeasible":
            # 当前前缀下取 0 不可行，故字典序最小解该位必为 1；
            # 上一次可行解必然已满足 x_i=1，直接沿用，无需再次求解。
            fixed[i] = 1
        else:
            return SolvedComponent(sites=comp_sites, decided=False)

    lex_haplotype = [0] + [fixed[i] for i in range(1, n)]

    # 第 3 步：可翻转性。每位点尝试在锚点为 0、代价最优的条件下取“字典序解的反值”。
    # 若可行，说明最优解集合中该位点既能取 0 也能取 1。
    flippable: List[int] = []
    for i in range(1, n):
        opposite = 1 - lex_haplotype[i]
        cm = build_component_model(
            comp_sites, comp_reads, fixed={i: opposite}, force_cost=opt_cost
        )
        kind, _ = run_solver(cm, time_limit, hint_x=lex_haplotype)
        if kind == "feasible":
            flippable.append(comp_sites[i])
        elif kind == "infeasible":
            continue
        else:
            return SolvedComponent(sites=comp_sites, decided=False)

    return SolvedComponent(
        sites=comp_sites,
        min_cost=opt_cost,
        haplotype=lex_haplotype,
        flippable_sites=sorted(flippable),
        decided=True,
    )


def phase_reads(
    n_sites: int, reads: List[Read], time_limit_seconds: Optional[float] = None
) -> List[SolvedComponent]:
    """整块求解入口：分块、逐块精确求解，不推断块间相位。"""
    time_limit = float(time_limit_seconds) if time_limit_seconds else DEFAULT_TIME_LIMIT
    return [
        solve_component(comp_sites, comp_reads, time_limit)
        for comp_sites, comp_reads in split_components(n_sites, reads)
    ]
