FROM sbubmi/segmentation_base:latest

# Clone this branch
RUN git clone -b subbranch_of_nscale_int https://github.com/SBU-BMI/pathomics_analysis.git && \
	cd pathomics_analysis && \
	cd nucleusSegmentation && \
	mkdir build && \
	cd build && \
	cmake ../src && \
	make -j4 && \
	cp app/main* /usr/local/bin/. && \
        cp ../script/mainAggregateFeatures.py  /usr/local/bin/.
