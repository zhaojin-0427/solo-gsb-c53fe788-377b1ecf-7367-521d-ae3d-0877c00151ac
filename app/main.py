"""FastAPI 入口：POST /phase。"""

from __future__ import annotations

from fastapi import FastAPI

from .phaser import UndecidedError, config_from_env, run_phasing
from .schemas import PhaseRequest, PhaseResponse

app = FastAPI(
    title="Diploid Short-Read Phasing Service",
    version="1.0.0",
    description=(
        "二倍体短读段定相服务：将读段归入两条互补单倍型，"
        "以不符观测的权重之和为代价，用 OR-Tools CP-SAT 求已证明的全局最小值。"
    ),
)


@app.get("/health")
def health() -> dict:
    return {"status": "ok"}


@app.post("/phase", response_model=PhaseResponse)
def phase(request: PhaseRequest) -> PhaseResponse:
    config = config_from_env()
    try:
        blocks = run_phasing(request, config)
    except UndecidedError as exc:
        # 未证明最优：明确返回 undecided，绝不给出伪最优结果
        return PhaseResponse(status="undecided", total_cost=None, blocks=None, detail=str(exc))
    return PhaseResponse(
        status="optimal",
        total_cost=sum(b.min_cost for b in blocks),
        blocks=blocks,
        detail=None,
    )
