FROM sbubmi/segmentation_base:latest

# Install analysis code
RUN git clone https://github.com/SBU-BMI/pathomics_analysis.git && \
	cd pathomics_analysis && \
	cd nucleusSegmentation && \
	mkdir build && \
	cd build && \
	cmake ../src && \
	make -j4 && \
	cp app/main* /usr/local/bin/. && \
        cp ../script/mainAggregateFeatures.py  /usr/local/bin/.
# Don't delete
	#cp ../script/mainAggregateFeatures.py  /usr/local/bin/. && \
	#cd /tmp && rm -rf pathomics_analysis*
