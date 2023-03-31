from __future__ import print_function

from dataclasses import dataclass

from typing import Dict,  List

from mzlib.attributes import AttributeManager, AttributeManagedProperty
from .utils import ensure_iter, flatten

SIMILAR_SPECTRUM_KEYS = "MS:1003263|similar spectrum keys"
SIMILAR_SPECTRUM_USI = "MS:1003264|similar spectrum USI"

CLUSTER_KEY = "MS:1003267|spectrum cluster key"

CLUSTER_MEMBERS_KEYS = "MS:1003268|spectrum cluster member spectrum keys"
CLUSTER_MEMBER_USI = "MS:1003269|spectrum cluster member USI"


@dataclass
class ClusterMemberRef:
    pass


@dataclass
class SpectrumRef(ClusterMemberRef):
    key: str


@dataclass
class USIRef(ClusterMemberRef):
    usi: str

    @property
    def key(self) -> str:
        return self.usi


class SpectrumCluster(AttributeManager):
    def __init__(self, attributes: List):
        super().__init__(attributes)

    key = AttributeManagedProperty[int](CLUSTER_KEY)
    _member_references = AttributeManagedProperty(CLUSTER_MEMBERS_KEYS)
    _cluster_member_usis = AttributeManagedProperty(CLUSTER_MEMBER_USI)

    @property
    def members(self) -> List[ClusterMemberRef]:
        internal_refs = [
            SpectrumRef(k) for k in flatten(ensure_iter(self._member_references))
        ]
        usi_members = [USIRef(k) for k in ensure_iter(self._cluster_member_usis)]
        return internal_refs + usi_members
