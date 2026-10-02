"""FastAPI 应用：二倍体短读段定相服务。"""
from __future__ import annotations

from fastapi import FastAPI

from .models import ComponentResult, PhaseRequest, PhaseResponse
from .solver import phase_reads

app = FastAPI(
    title="Diploid Short-Read Phasing Service",
    description="基于 OR-Tools CP-SAT 的加权最小纠错（MEC）精确定相服务",
    version="1.0.0",
)


@app.get("/health")
def health() -> dict:
    return {"status": "ok"}


@app.post("/phase", response_model=PhaseResponse)
def phase(request: PhaseRequest) -> PhaseResponse:
    solved = phase_reads(
        n_sites=len(request.sites),
        reads=request.reads,
        time_limit_seconds=request.time_limit_seconds,
    )

    all_decided = all(component.decided for component in solved)
    components = [
        ComponentResult(
            status="optimal" if component.decided else "undecided",
            sites=component.sites,
            min_cost=component.min_cost if component.decided else None,
            haplotype=component.haplotype if component.decided else None,
            flippable_sites=component.flippable_sites if component.decided else None,
        )
        for component in solved
    ]
    return PhaseResponse(status="optimal" if all_decided else "undecided", components=components)
