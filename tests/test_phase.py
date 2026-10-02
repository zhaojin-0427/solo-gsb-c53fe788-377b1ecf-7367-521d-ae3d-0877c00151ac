"""接口与算法测试。

正确性测试用暴力枚举作为参照：对每个随机小实例，枚举块内所有
单倍型（首位固定 0），独立计算最小代价、字典序最小单倍型和可翻转
位点，与服务返回结果比对。
"""

from __future__ import annotations

import itertools
import random

import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


# ---------------------------------------------------------------------------
# 工具
# ---------------------------------------------------------------------------

def make_payload(sites, reads):
    return {
        "sites": [{"id": s[0], "position": s[1]} for s in sites],
        "reads": [
            {"observations": [{"site": o[0], "allele": o[1], "weight": o[2]} for o in obs]}
            for obs in reads
        ],
    }


def brute_force_block(site_ids, block_reads):
    """对单个块暴力求解。block_reads: [[(site_id, allele, weight), ...], ...]"""
    n = len(site_ids)
    idx = {s: i for i, s in enumerate(site_ids)}
    best_cost, best_haps = None, []
    for bits in itertools.product([0, 1], repeat=n - 1):
        hap = (0,) + bits  # 首位固定 0
        comp = tuple(1 - b for b in hap)
        cost = 0
        for obs in block_reads:
            e0 = sum(w for s, a, w in obs if hap[idx[s]] != a)
            e1 = sum(w for s, a, w in obs if comp[idx[s]] != a)
            cost += min(e0, e1)
        if best_cost is None or cost < best_cost:
            best_cost, best_haps = cost, [hap]
        elif cost == best_cost:
            best_haps.append(hap)
    lex = min(best_haps)
    flippable = sorted(
        {i for hap in best_haps for i in range(n) if hap[i] != lex[i]}
    )
    return best_cost, "".join(map(str, lex)), [site_ids[i] for i in flippable]


# ---------------------------------------------------------------------------
# 422 校验
# ---------------------------------------------------------------------------

def test_duplicate_site_id_422():
    payload = make_payload([("s1", 100), ("s1", 200)], [])
    assert client.post("/phase", json=payload).status_code == 422


def test_unsorted_positions_422():
    payload = make_payload([("s1", 200), ("s2", 100)], [])
    assert client.post("/phase", json=payload).status_code == 422


def test_duplicate_position_422():
    payload = make_payload([("s1", 100), ("s2", 100)], [])
    assert client.post("/phase", json=payload).status_code == 422


def test_read_with_fewer_than_two_observations_422():
    payload = make_payload([("s1", 100), ("s2", 200)], [[("s1", 0, 30)]])
    assert client.post("/phase", json=payload).status_code == 422


def test_read_with_unknown_site_422():
    payload = make_payload(
        [("s1", 100), ("s2", 200)], [[("s1", 0, 30), ("sX", 1, 30)]]
    )
    assert client.post("/phase", json=payload).status_code == 422


def test_read_covering_same_site_twice_422():
    payload = make_payload(
        [("s1", 100), ("s2", 200)], [[("s1", 0, 30), ("s1", 1, 30)]]
    )
    assert client.post("/phase", json=payload).status_code == 422


@pytest.mark.parametrize("allele", [-1, 2, 1.5, True, "0"])
def test_invalid_allele_422(allele):
    payload = make_payload(
        [("s1", 100), ("s2", 200)], [[("s1", allele, 30), ("s2", 1, 30)]]
    )
    assert client.post("/phase", json=payload).status_code == 422


@pytest.mark.parametrize("weight", [0, 61, -5, 2.5, True, "30"])
def test_invalid_weight_422(weight):
    payload = make_payload(
        [("s1", 100), ("s2", 200)], [[("s1", 0, weight), ("s2", 1, 30)]]
    )
    assert client.post("/phase", json=payload).status_code == 422


def test_empty_sites_422():
    assert client.post("/phase", json={"sites": [], "reads": []}).status_code == 422


# ---------------------------------------------------------------------------
# 正确性
# ---------------------------------------------------------------------------

def test_two_blocks_with_flippable():
    """README 中的示例：两个块，第二块存在等权冲突导致 s5 可翻转。"""
    payload = make_payload(
        [("s1", 100), ("s2", 200), ("s3", 300), ("s4", 400), ("s5", 500)],
        [
            [("s1", 0, 30), ("s2", 0, 30)],
            [("s1", 0, 10), ("s2", 1, 10)],
            [("s2", 0, 30), ("s3", 1, 30)],
            [("s4", 1, 20), ("s5", 1, 20)],
            [("s4", 1, 20), ("s5", 0, 20)],
        ],
    )
    resp = client.post("/phase", json=payload)
    assert resp.status_code == 200
    body = resp.json()
    assert body["status"] == "optimal"
    assert body["total_cost"] == 30
    assert body["blocks"] == [
        {
            "index": 0,
            "site_ids": ["s1", "s2", "s3"],
            "min_cost": 10,
            "haplotype": "001",
            "flippable_site_ids": [],
        },
        {
            "index": 1,
            "site_ids": ["s4", "s5"],
            "min_cost": 20,
            "haplotype": "00",
            "flippable_site_ids": ["s5"],
        },
    ]


def test_isolated_site_is_singleton_block():
    payload = make_payload(
        [("s1", 100), ("s2", 200), ("s3", 300)],
        [[("s1", 1, 30), ("s2", 1, 30)]],
    )
    body = client.post("/phase", json=payload).json()
    assert body["status"] == "optimal"
    assert body["blocks"][0] == {
        "index": 0,
        "site_ids": ["s1", "s2"],
        "min_cost": 0,
        "haplotype": "00",
        "flippable_site_ids": [],
    }
    assert body["blocks"][1] == {
        "index": 1,
        "site_ids": ["s3"],
        "min_cost": 0,
        "haplotype": "0",
        "flippable_site_ids": [],
    }


def test_empty_reads_all_singletons():
    payload = make_payload([("a", 1), ("b", 2)], [])
    body = client.post("/phase", json=payload).json()
    assert body["status"] == "optimal"
    assert body["total_cost"] == 0
    assert [b["haplotype"] for b in body["blocks"]] == ["0", "0"]


def test_interleaved_blocks():
    """读段交错覆盖时仍按共覆盖正确分块：{s1,s3} 与 {s2,s4}。"""
    payload = make_payload(
        [("s1", 100), ("s2", 200), ("s3", 300), ("s4", 400)],
        [
            [("s1", 0, 30), ("s3", 0, 30)],
            [("s2", 1, 30), ("s4", 1, 30)],
        ],
    )
    body = client.post("/phase", json=payload).json()
    assert body["status"] == "optimal"
    assert [b["site_ids"] for b in body["blocks"]] == [["s1", "s3"], ["s2", "s4"]]
    assert [b["haplotype"] for b in body["blocks"]] == ["00", "00"]


def test_complement_symmetry_broken():
    """单条读段 (1,1)：互补解 (0,0)/(1,1) 等优，首位固定 0 后应返回 00。"""
    payload = make_payload(
        [("s1", 100), ("s2", 200)], [[("s1", 1, 30), ("s2", 1, 30)]]
    )
    body = client.post("/phase", json=payload).json()
    assert body["blocks"][0]["haplotype"] == "00"
    assert body["blocks"][0]["min_cost"] == 0


def test_random_small_instances_against_brute_force():
    rng = random.Random(20261002)
    for trial in range(60):
        n_sites = rng.randint(2, 7)
        sites = [(f"s{i}", (i + 1) * 100) for i in range(n_sites)]
        reads = []
        for _ in range(rng.randint(1, 10)):
            k = rng.randint(2, min(4, n_sites))
            chosen = rng.sample([s[0] for s in sites], k)
            reads.append(
                [(c, rng.randint(0, 1), rng.randint(1, 60)) for c in chosen]
            )
        resp = client.post("/phase", json=make_payload(sites, reads))
        assert resp.status_code == 200
        body = resp.json()
        assert body["status"] == "optimal", f"trial {trial}: {body}"
        # 按返回的块分组，逐块与暴力解比对
        read_by_block = {b["index"]: [] for b in body["blocks"]}
        block_of = {}
        for b in body["blocks"]:
            for sid in b["site_ids"]:
                block_of[sid] = b["index"]
        for obs in reads:
            read_by_block[block_of[obs[0][0]]].append(obs)
        for b in body["blocks"]:
            cost, hap, flip = brute_force_block(
                b["site_ids"], read_by_block[b["index"]]
            )
            assert b["min_cost"] == cost, f"trial {trial} block {b}"
            assert b["haplotype"] == hap, f"trial {trial} block {b}"
            assert b["flippable_site_ids"] == flip, f"trial {trial} block {b}"
        assert body["total_cost"] == sum(b["min_cost"] for b in body["blocks"])


# ---------------------------------------------------------------------------
# undecided 行为
# ---------------------------------------------------------------------------

def test_undecided_when_optimality_not_proven(monkeypatch):
    """模拟求解器无法证明最优时，必须返回 undecided 且不给出块结果。"""
    import app.main as main_module

    def fake_run_phasing(request, config):
        raise main_module.UndecidedError("primary solve did not prove optimality (status=FEASIBLE)")

    monkeypatch.setattr(main_module, "run_phasing", fake_run_phasing)
    payload = make_payload(
        [("s1", 100), ("s2", 200)], [[("s1", 0, 30), ("s2", 1, 30)]]
    )
    resp = client.post("/phase", json=payload)
    assert resp.status_code == 200
    body = resp.json()
    assert body["status"] == "undecided"
    assert body["blocks"] is None
    assert body["total_cost"] is None
    assert "did not prove optimality" in body["detail"]


def test_health():
    assert client.get("/health").json() == {"status": "ok"}
