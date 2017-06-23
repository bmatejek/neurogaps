#
# Makefile for GAPS apps
#

opt:
	$(MAKE) target "TARGET=$@"

_debug:
	$(MAKE) target "TARGET=$@"

clean:
	$(MAKE) target "TARGET=$@"

release:
	mkdir -p ../release/apps
	cp Makefile ../release/apps
	$(MAKE) target "TARGET=$@"

target: 
	cd aff2prob; $(MAKE) $(TARGET)
	cd affinitymap; $(MAKE) $(TARGET)
	cd augmentdata; $(MAKE) $(TARGET)
	cd averagecommute; $(MAKE) $(TARGET)
	cd boundarygrowth; $(MAKE) $(TARGET)
	cd boundarypaths; $(MAKE) $(TARGET)
	cd boundarymerge; $(MAKE) $(TARGET)
	cd createdata; $(MAKE) $(TARGET)
	cd createforest; $(MAKE) $(TARGET)
	cd deflateraw; $(MAKE) $(TARGET)
	cd derivative; $(MAKE) $(TARGET)
	cd dijkstra; $(MAKE) $(TARGET)
	cd dijkstraedge; $(MAKE) $(TARGET)
	cd dijkstrapost; $(MAKE) $(TARGET)
	cd dijkstrastat; $(MAKE) $(TARGET)
	cd dijkstraview; $(MAKE) $(TARGET)
	cd entropy; $(MAKE) $(TARGET)
	cd evaluatedist; $(MAKE) $(TARGET)
	cd facilityloc; $(MAKE) $(TARGET)
	cd forestresults; $(MAKE) $(TARGET)
	cd genfeatures; $(MAKE) $(TARGET)
	cd globaldata; $(MAKE) $(TARGET)
	cd graphcut; $(MAKE) $(TARGET)
	cd graphcutbinary; $(MAKE) $(TARGET)
	cd graphcutstat; $(MAKE) $(TARGET)
	cd graphcutunary; $(MAKE) $(TARGET)
	cd graphcutview; $(MAKE) $(TARGET)
	cd grd2meta; $(MAKE) $(TARGET)
	cd grdview; $(MAKE) $(TARGET)
	cd hierarchical; $(MAKE) $(TARGET)
	cd hierordering; $(MAKE) $(TARGET)
	cd kmedian; $(MAKE) $(TARGET)
	cd lpboundary; $(MAKE) $(TARGET)
	cd meta2grd; $(MAKE) $(TARGET)
	cd metaimage; $(MAKE) $(TARGET)
	cd metamerge; $(MAKE) $(TARGET)
	cd metaslice; $(MAKE) $(TARGET)
	cd metasplit; $(MAKE) $(TARGET)
	cd ml2pred; $(MAKE) $(TARGET)
	cd montecarlo; $(MAKE) $(TARGET)
	cd neurondata; $(MAKE) $(TARGET)
	cd neuronfilter; $(MAKE) $(TARGET)
	cd neuronread; $(MAKE) $(TARGET)
	cd neuronstat; $(MAKE) $(TARGET)
	cd neurontest; $(MAKE) $(TARGET)
	cd neurontruth; $(MAKE) $(TARGET)
	cd neuronview; $(MAKE) $(TARGET)
	cd outputview; $(MAKE) $(TARGET)
	cd ordermerge; $(MAKE) $(TARGET)
	cd pathvisual; $(MAKE) $(TARGET)
	cd postcreatedata; $(MAKE) $(TARGET)
	cd postcreateforest; $(MAKE) $(TARGET)
	cd postrandomforest; $(MAKE) $(TARGET)
	cd postprocess; $(MAKE) $(TARGET)
	cd randomforest; $(MAKE) $(TARGET)
	cd randomhier; $(MAKE) $(TARGET)
	cd refinemeta; $(MAKE) $(TARGET)
	cd rfview; $(MAKE) $(TARGET)
	cd setcompare; $(MAKE) $(TARGET)
	cd skel2grd; $(MAKE) $(TARGET)
	cd skelcheck; $(MAKE) $(TARGET)
	cd skelstat; $(MAKE) $(TARGET)
	cd statistics; $(MAKE) $(TARGET)
	cd subsection; $(MAKE) $(TARGET)
	cd super2graph $(MAKE) $(TARGET)
	cd testmetric; $(MAKE) $(TARGET)
	cd txt2neuron; $(MAKE) $(TARGET)
	cd vesselstat; $(MAKE) $(TARGET)
	cd vesselview; $(MAKE) $(TARGET)
	cd voxelpair; $(MAKE) $(TARGET)
