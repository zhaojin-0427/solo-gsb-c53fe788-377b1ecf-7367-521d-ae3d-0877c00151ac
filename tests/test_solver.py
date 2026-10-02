"""定相逻辑测试：连通块、加权 MEC、字典序最小解、可翻转位点，全部与暴力枚举对拍。"""
from __future__ import annotations

import itertools
import random
from typing import List, Tuple

import pytest

from app.models import Observation, PhaseRequest, Read
from app.solver import phase_reads, solve_component, split_components


def make_read(rid: str, obs: List[Tuple[int, int, int]]) -> Read:
    return Read(
        id=rid,
        observations=[Observation(site=s, allele=a, weight=w) for s, a, w in obs],
    )


def brute_force(n_sites: int, reads: List[Read]):
    """枚举 (x, y)，返回 (最优代价, 锚定 x0=0 的全部最优 x 列表)。"""""
    comps = split_components(n_sites, reads)
    per_component = []
    for comp_sites, comp_reads in comps:
        m = len(comp_sites)
        best = None
        best_xs = []
        for x_bits in itertools.product([0, 1], repeat=m):
            if x_bits[0] != 0:
                continue
            for y_bits in itertools.product([0, 1], repeat=len(comp_reads)):
                cost = 0
                for r_idx, read in enumerate(comp_reads):
                    li = {g: i for i, g in enumerate(comp_sites)}
                    for obs in read.observations:
                        mismatch = x_bits[li[obs.site]] ^ obs.allele ^ y_bits[r_idx]
                        cost += mismatch * obs.weight
                if best is None or cost < best:
                    best = cost
                    best_xs = [x_bits]
                elif cost == best:
                    best_xs.append(x_bits)
        unique_xs = sorted(set(best_xs))
        flippable_local = [i for i in range(m) if len({x[i] for x in unique_xs}) == 2]
        per_component.append(
            (comp_sites, best, unique_xs[0] if unique_xs else None,
             [comp_sites[i] for i in flippable_local])
        )
    return per_component


def assert_matches_bruteforce(n_sites: int, reads: List[Read], time_limit: float = 10.0):
    solved = phase_reads(n_sites, reads, time_limit)
    expected = brute_force(n_sites, reads)
    assert len(solved) == len(expected)
    for got, (sites, cost, lex_hap, flip) in zip(solved, expected):
        assert got.decided is True
        assert got.sites == sites
        assert got.min_cost == cost
        assert got.haplotype == list(lex_hap)
        assert got.flippable_sites == flip


def test_simple_two_site():
    # r1 同时观察到 0、1：x0=0, x1=1 时零代价。
    reads = [make_read("r1", [(0, 0, 1), (1, 1, 1)])]
    assert_matches_bruteforce(2, reads)
    solved = phase_reads(2, reads, 10.0)
    assert solved[0].haplotype == [0, 1]
    assert solved[0].min_cost == 0
    assert solved[0].flippable_sites == []


def test_weighted_conflict():
    # 两条读段互相矛盾：r1 说 00，r2 说 01，代价小者输。
    reads = [
        make_read("r1", [(0, 0, 1), (1, 0, 1)]),
        make_read("r2", [(0, 0, 1), (1, 1, 1)]),
    ]
    solved = phase_reads(2, reads, 10.0)
    assert solved[0].min_cost == 1
    assert solved[0].haplotype == [0, 0]  # 字典序：翻转后 x1=1 等价，取 0
    # x1=0 与 x1=1 都是最优 -> 可翻转
    assert solved[0].flippable_sites == [1]


def test_flippable_pair():
    # 单条读段只覆盖位点 1、2（位点 0 孤立），块内首位（位点 1）锚定为 0，
    # 位点 2 被读段约束为 1，无可翻转；再构造一个完全自由翻转的块。
    reads = [make_read("r1", [(2, 0, 1), (3, 1, 1)])]
    solved = phase_reads(4, reads, 10.0)
    # 三个块：{0} 孤立, {1} 孤立, {2,3}
    assert [c.sites for c in solved] == [[0], [1], [2, 3]]
    assert solved[0].min_cost == 0 and solved[0].haplotype == [0]
    assert solved[2].haplotype == [0, 1]
    assert solved[2].flippable_sites == []


def test_free_block_all_flippable():
    # 三条位点成块但读段信息自洽且整体翻转等价：构造两对“成对翻转”自由度。
    # r1: 0-1 为 00；r2: 1-2 为 11（相对 r1 在 1 上取反），
    # 这样最优唯一 [0,0,1]。再构造真正的自由块：仅有一条读段覆盖 0,1 说 01，
    # 另一条读段覆盖 2,3 说 01，两个块，各自无翻转。
    reads = [
        make_read("a", [(0, 0, 1), (1, 1, 1)]),
        make_read("b", [(2, 0, 1), (3, 1, 1)]),
    ]
    assert_matches_bruteforce(4, reads)


def test_paired_flip_is_flippable():
    # 构造位点 1、2 必须成对翻转的块（锚点位点 0 固定为 0）：
    #   r1 覆盖 (0,1) 观测 01：x1 必须与 x0 不同，即 x1=1；
    #   r2 覆盖 (0,2) 观测 01：x2 必须与 x0 不同，即 x2=1；
    #   r3 覆盖 (1,2) 观测 00（权重 1）：要求 x1==x2，零代价时成立。
    # 另有弱权重读段 r4 覆盖 (0,1) 观测 00（权重 1）：与 r1 矛盾，
    # 违背 r1/r2 或违背 r3+r4 代价均为 2，于是最优解集合为
    #   [0,1,1]（代价 1: 违背 r4）与 [0,0,0]（代价 2？）...
    # 为保证严格成对等代价，改用等权重冲突：
    #   r1 (0,1)=01 w=1, r2 (0,2)=01 w=1, r3 (1,2)=00 w=1,
    #   r4 (0,1)=00 w=1, r5 (0,2)=00 w=1。
    # 解 [0,1,1]：r4、r2 与 r5? r5 要求 x2=x0 -> 违背；代价 r4+r5=2。
    # 解 [0,0,0]：r1、r2 违背，r3/r4/r5 满足；代价 2。
    # 故位点 1、2 均可翻转（且成对），字典序最小解为 [0,0,0]。
    reads = [
        make_read("r1", [(0, 0, 1), (1, 1, 1)]),
        make_read("r2", [(0, 0, 1), (2, 1, 1)]),
        make_read("r3", [(1, 0, 1), (2, 0, 1)]),
        make_read("r4", [(0, 0, 1), (1, 0, 1)]),
        make_read("r5", [(0, 0, 1), (2, 0, 1)]),
    ]
    solved = phase_reads(3, reads, 10.0)
    assert solved[0].min_cost == 2
    assert solved[0].haplotype == [0, 0, 0]
    assert solved[0].flippable_sites == [1, 2]
    assert_matches_bruteforce(3, reads)


def test_weights_break_ties_and_boundary():
    # 权重 1 与 60：60 的观测决不能被违背，除非必要。
    reads = [
        make_read("big", [(0, 0, 60), (1, 1, 60)]),
        make_read("small", [(0, 0, 1), (1, 0, 1)]),
    ]
    solved = phase_reads(2, reads, 10.0)
    assert solved[0].min_cost == 1
    assert solved[0].haplotype == [0, 1]
    assert solved[0].flippable_sites == []


def test_random_against_bruteforce():
    rng = random.Random(20261002)
    for trial in range(40):
        n = rng.randint(1, 5)
        n_reads = rng.randint(1, 5)
        reads = []
        for rid in range(n_reads):
            k = rng.randint(2, min(4, n)) if n >= 2 else 2
            if n < 2:
                break
            sites = rng.sample(range(n), k)
            obs = [(s, rng.randint(0, 1), rng.randint(1, 60)) for s in sites]
            reads.append(make_read(f"r{rid}", obs))
        if n < 2:
            continue
        assert_matches_bruteforce(n, reads)


# ---------------------- HTTP 层校验 ---------------------- #
def test_validation_unsorted_sites():
    req = PhaseRequest.model_validate  # 仅为可读性
    from fastapi.testclient import TestClient
    from app.main import app

    client = TestClient(app)
    body = {
        "sites": [{"position": 10}, {"position": 5}],
        "reads": [{"id": "r1", "observations": [
            {"site": 0, "allele": 0, "weight": 1},
            {"site": 1, "allele": 1, "weight": 1},
        ]}],
    }
    resp = client.post("/phase", json=body)
    assert resp.status_code == 422


def test_validation_duplicate_positions():
    from fastapi.testclient import TestClient
    from app.main import app

    client = TestClient(app)
    body = {
        "sites": [{"position": 7}, {"position": 7}],
        "reads": [{"id": "r1", "observations": [
            {"site": 0, "allele": 0, "weight": 1},
            {"site": 1, "allele": 1, "weight": 1},
        ]}],
    }
    resp = client.post("/phase", json=body)
    assert resp.status_code == 422


def test_validation_illegal_reads():
    from fastapi.testclient import TestClient
    from app.main import app

    client = TestClient(app)
    base_sites = [{"position": 0}, {"position": 1}, {"position": 2}]

    # 读段只覆盖 1 个位点
    body = {"sites": base_sites, "reads": [{"id": "r1", "observations": [
        {"site": 0, "allele": 0, "weight": 1}]}]}
    assert client.post("/phase", json=body).status_code == 422

    # 非法等位
    body = {"sites": base_sites, "reads": [{"id": "r1", "observations": [
        {"site": 0, "allele": 2, "weight": 1},
        {"site": 1, "allele": 0, "weight": 1}]}]}
    assert client.post("/phase", json=body).status_code == 422

    # 权重越界（0 与 61）
    for bad_weight in (0, 61):
        body = {"sites": base_sites, "reads": [{"id": "r1", "observations": [
            {"site": 0, "allele": 0, "weight": bad_weight},
            {"site": 1, "allele": 0, "weight": 1}]}]}
        assert client.post("/phase", json=body).status_code == 422

    # site 下标越界
    body = {"sites": base_sites, "reads": [{"id": "r1", "observations": [
        {"site": 0, "allele": 0, "weight": 1},
        {"site": 5, "allele": 0, "weight": 1}]}]}
    assert client.post("/phase", json=body).status_code == 422

    # 同读段重复位点
    body = {"sites": base_sites, "reads": [{"id": "r1", "observations": [
        {"site": 0, "allele": 0, "weight": 1},
        {"site": 0, "allele": 1, "weight": 1}]}]}
    assert client.post("/phase", json=body).status_code == 422


def test_end_to_end_api():
    from fastapi.testclient import TestClient
    from app.main import app

    client = TestClient(app)
    assert client.get("/health").json() == {"status": "ok"}

    body = {
        "sites": [{"position": 100}, {"position": 200}, {"position": 300}],
        "reads": [
            {"id": "r1", "observations": [
                {"site": 0, "allele": 0, "weight": 30},
                {"site": 1, "allele": 1, "weight": 30}]},
            {"id": "r2", "observations": [
                {"site": 1, "allele": 1, "weight": 20},
                {"site": 2, "allele": 0, "weight": 20}]},
        ],
    }
    resp = client.post("/phase", json=body)
    assert resp.status_code == 200
    data = resp.json()
    assert data["status"] == "optimal"
    assert len(data["components"]) == 1
    comp = data["components"][0]
    assert comp["sites"] == [0, 1, 2]
    assert comp["min_cost"] == 0
    assert comp["haplotype"] == [0, 1, 0]
    assert comp["flippable_sites"] == []


def test_multi_component_response():
    from fastapi.testclient import TestClient
    from app.main import app

    client = TestClient(app)
    body = {
        "sites": [{"position": p} for p in (10, 20, 30, 40)],
        "reads": [
            {"id": "a", "observations": [
                {"site": 0, "allele": 0, "weight": 1},
                {"site": 1, "allele": 1, "weight": 1}]},
            {"id": "b", "observations": [
                {"site": 2, "allele": 1, "weight": 1},
                {"site": 3, "allele": 0, "weight": 1}]},
        ],
    }
    data = client.post("/phase", json=body).json()
    assert data["status"] == "optimal"
    assert [c["sites"] for c in data["components"]] == [[0, 1], [2, 3]]
    assert data["components"][0]["haplotype"] == [0, 1]
    assert data["components"][1]["haplotype"] == [0, 1]
