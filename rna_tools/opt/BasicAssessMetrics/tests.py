#!/usr/bin/env python
# -*- coding: utf-8 -*-

# calculate InteractionNetworkFidelity and Deformation Index for RNA structures
# need to have MA-annotate in the directory or set in mcannotate.py
from BasicAssessMetrics import normalize_structure, \
                               interaction_network_fidelity, \
                               calc_RMSD

## fn = "example/14_BujnickiPreExp_2.pdb"
## i1 = "example/14_BujnickiPreExp_2.index"
## i2 = i1
## rmsd, DI_ALL, INF_ALL, INF_WC, INF_NWC,INF_STACK = InteractionNetworkFidelity(fn, i1, fn, i2)

## print '14_ChenPostExp_2, rmsd', rmsd
## print "  DI_ALL:", DI_ALL
## print "  INF_ALL:", INF_ALL


fn = "example/14_BujnickiPreExp_2.pdb"
i1 = ""
i2 = ""
rmsd, DI_ALL, INF_ALL, INF_WC, INF_NWC, INF_STACK =  interaction_network_fidelity(fn, None, fn, None)

print "  rmsd:", rmsd
print "  DI_ALL:", DI_ALL
print "  INF_ALL:", INF_ALL
