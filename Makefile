CXX = clang++ -Wall -Wextra -O3 --std=c++11 

all: volume2positionList label2color


volume2positionList: main.cpp tomogram.hpp bilateral.hpp gauss_filter.hpp common.hpp homogenator.hpp cSegmentationByOtsu.h watershed_segmentation.hpp colors.hpp
	$(CXX) main.cpp -o volume2positionList

label2color: label2color.cpp tomogram.hpp bilateral.hpp gauss_filter.hpp common.hpp homogenator.hpp cSegmentationByOtsu.h watershed_segmentation.hpp colors.hpp
	$(CXX) label2color.cpp -o label2color

clean:
	rm volume2positionList
	rm label2color
