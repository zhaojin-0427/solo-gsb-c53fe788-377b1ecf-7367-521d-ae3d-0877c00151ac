"""POST /phase 的请求与响应数据模型。

输入与输出的等位取值均为 0/1；权重是 1~60 的整数。
所有字段校验失败都由 FastAPI 转成 HTTP 422 响应。
"""
from __future__ import annotations

from typing import List, Optional

from pydantic import BaseModel, ConfigDict, Field, model_validator


class PhaseSite(BaseModel):
    """杂合双等位位点，只需给出基因组位置。"""

    model_config = ConfigDict(strict=True, extra="forbid")

    position: int = Field(..., ge=0, description="基因组位置（非负整数），请求中必须唯一且严格递增")


class Observation(BaseModel):
    """一条读段在某个已知位点上的等位观测。"""

    model_config = ConfigDict(strict=True, extra="forbid")

    site: int = Field(..., ge=0, description="观测所在位点在请求 sites 列表中的下标（从 0 开始）")
    allele: int = Field(..., ge=0, le=1, description="观测到的等位：0 或 1")
    weight: int = Field(..., ge=1, le=60, description="该观测的整数权重（1~60）")


class Read(BaseModel):
    """一条覆盖至少两个已知位点的读段。"""

    model_config = ConfigDict(strict=True, extra="forbid")

    id: str = Field(..., min_length=1, max_length=128)
    observations: List[Observation] = Field(..., min_length=2)


class PhaseRequest(BaseModel):
    """POST /phase 请求体。"""

    model_config = ConfigDict(strict=True, extra="forbid")

    sites: List[PhaseSite] = Field(..., min_length=1)
    reads: List[Read] = Field(..., min_length=1)
    time_limit_seconds: Optional[int] = Field(
        default=None,
        ge=1,
        le=120,
        description="每次求解器调用的最长秒数（1~120）；默认 20 秒",
    )

    @model_validator(mode="after")
    def _validate_request(self) -> "PhaseRequest":
        # 位点：按位置严格递增，等价于“已排序 + 唯一”。
        positions = [s.position for s in self.sites]
        if any(positions[i] >= positions[i + 1] for i in range(len(positions) - 1)):
            raise ValueError("sites 必须按 position 严格递增且不得包含重复位置")

        n = len(self.sites)
        seen_read_ids: set[str] = set()
        for read in self.reads:
            if read.id in seen_read_ids:
                raise ValueError(f"读段 id 重复：{read.id!r}")
            seen_read_ids.add(read.id)

            sites_in_read = [obs.site for obs in read.observations]
            if any(site >= n for site in sites_in_read):
                raise ValueError(f"读段 {read.id!r} 的 observation.site 越界（位点总数为 {n}）")
            if len(set(sites_in_read)) != len(sites_in_read):
                raise ValueError(f"读段 {read.id!r} 在同一位点上存在重复观测")
        return self


class ComponentResult(BaseModel):
    """单个共覆盖连通块的结果。

    status 为 "undecided" 时只给出块内位点下标和状态，不给出任何未被证明的数值。
    """

    status: str = Field(..., description="optimal 或 undecided")
    sites: List[int] = Field(..., description="块内位点的全局下标（按位置排序）")
    min_cost: Optional[int] = Field(default=None, description="该块最优代价（不符观测权重之和）")
    haplotype: Optional[List[int]] = Field(
        default=None, description="按位点顺序给出的字典序最小最优单倍型（0/1）"
    )
    flippable_sites: Optional[List[int]] = Field(
        default=None,
        description=(
            "在所有最优解中可以翻转的位点（全局下标）。"
            "即同时存在最优解使该位点取 0 与取 1；不要求其余位点保持不变"
        ),
    )


class PhaseResponse(BaseModel):
    status: str = Field(..., description="optimal：全部连通块已证明最优；undecided：存在未决块")
    components: List[ComponentResult]
