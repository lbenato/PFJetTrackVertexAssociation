# PFJetTrackVertexAssociation

```bash
mkdir PF_track_vertex_association
cmsrel CMSSW_9_4_1
cd CMSSW_9_4_1/src
cmsenv
git cms-init

#clone repo
mk Analyzer
cd Analyzer
git clone https://github.com/lbenato/PFJetTrackVertexAssociation.git
cd PFJetTrackVertexAssociation
scram b -j 32
```
