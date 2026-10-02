"""请求 / 响应的数据模型与校验规则。

所有校验失败（重复位点、未排序、非法读段等）都会抛出 ValueError，
由 FastAPI 统一转换为 HTTP 422 响应。
"""

from __future__ import annotations

from typing import List, Optional

from pydantic import BaseModel, Field, model_validator


class Site(BaseModel):
    """一个杂合双等位位点。"""

    id: str = Field(min_length=1, description="位点标识，请求内唯一")
    position: int = Field(description="基因组坐标，请求中必须严格递增")


class Observation(BaseModel):
    """读段在某个位点上的一次观测。"""

    site: str = Field(min_length=1, description="位点 id，必须出现在 sites 中")
    allele: int = Field(ge=0, le=1, strict=True, description="等位取值，只能为 0 或 1")
    weight: int = Field(ge=1, le=60, strict=True, description="观测权重，1~60 的整数")


class Read(BaseModel):
    """一条短读段，至少覆盖两个不同的已知位点。"""

    id: Optional[str] = Field(default=None, description="可选的读段标识，仅用于客户端追踪")
    observations: List[Observation] = Field(
        min_length=2, description="至少两条观测，且位点互不重复"
    )


class PhaseRequest(BaseModel):
    sites: List[Site] = Field(min_length=1, description="按位置严格递增排列的唯一杂合位点")
    reads: List[Read] = Field(default_factory=list, description="读段列表，可以为空")

    @model_validator(mode="after")
    def _validate(self) -> "PhaseRequest":
        # 位点 id 唯一
        ids = [s.id for s in self.sites]
        if len(set(ids)) != len(ids):
            dup = sorted({i for i in ids if ids.count(i) > 1})
            raise ValueError(f"duplicate site id(s): {dup}")
        # 位置严格递增（既保证唯一，也保证已按位置排序）
        positions = [s.position for s in self.sites]
        if any(b <= a for a, b in zip(positions, positions[1:])):
            raise ValueError("site positions must be strictly increasing (sorted and unique)")

        known = set(ids)
        for r_idx, read in enumerate(self.reads):
            seen = set()
            for obs in read.observations:
                if obs.site not in known:
                    raise ValueError(
                        f"reads[{r_idx}] references unknown site {obs.site!r}"
                    )
                if obs.site in seen:
                    raise ValueError(
                        f"reads[{r_idx}] covers site {obs.site!r} more than once"
                    )
                seen.add(obs.site)
        return self


class BlockResult(BaseModel):
    """一个连通块的定相结果。"""

    index: int = Field(description="块序号，按块内最小位点位置排序")
    site_ids: List[str] = Field(description="块内位点 id，与请求中的位点顺序一致")
    min_cost: int = Field(description="该块的最小代价（不符观测的权重之和）")
    haplotype: str = Field(
        description="字典序最小的最优单倍型（0/1 字符串，与 site_ids 对齐；首位固定为 0）"
    )
    flippable_site_ids: List[str] = Field(
        description="在所有最优解中可取另一等位的位点（块首位因对称性固定，永不包含）"
    )


class PhaseResponse(BaseModel):
    status: str = Field(description='"optimal"：已证明全局最优；"undecided"：未能在时限内证明最优')
    total_cost: Optional[int] = Field(default=None, description="各块最小代价之和")
    blocks: Optional[List[BlockResult]] = Field(
        default=None, description="各连通块结果；status 为 undecided 时为 null"
    )
    detail: Optional[str] = Field(default=None, description="undecided 时的原因说明")
