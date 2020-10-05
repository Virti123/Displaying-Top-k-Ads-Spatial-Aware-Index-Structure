//C++ Implementation of Top-k geometric intersection queries motivated by ad auctions.
#include <chrono> 
#include <fstream>
#include<iostream>
#include <boost/config.hpp>
#include <boost/version.hpp>
#include <fstream>
#include <random>
#include <math.h>
#include <unistd.h>

// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_hierarchy_2.h>
#include <CGAL/Apollonius_graph_filtered_traits_2.h>
#include <CGAL/point_generators_2.h>
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
//#include <CGAL/IO/WKT.h>
#endif

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>

// Qt headers
#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>

#include <CGAL/Qt/ApolloniusGraphGraphicsItem.h>
#include<CGAL/Qt/GraphicsViewCircleInput.h>
#include<CGAL/Voronoi_diagram_2/Apollonius_graph_nearest_site_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Qt/Converter.h>

// GraphicsView items and event filters (input classes)
//#include <CGAL/Qt/ApolloniusGraphGraphicsItem.h>
//#include <CGAL/Qt/GraphicsViewCircleInput.h>

// for viewportsBbox
#include <CGAL/Qt/utility.h>
#include<CGAL/Apollonius_graph_adaptation_traits_2.h>
// the two base classes
//#include "ui_Apollonius_graph_2.h"
#include <CGAL/Qt/DemosMainWindow.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Random.h>
#include <CGAL/Apollonius_graph_adaptation_policies_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <QPointF>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>



typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_2 Point_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;

typedef CGAL::Apollonius_graph_filtered_traits_2<K,CGAL::Integral_domain_without_division_tag>  Gt;

typedef Gt::Point_2                           Point_2;
typedef K::Circle_2                           Circle_2;
typedef Gt::Site_2                            Apollonius_site_2;
typedef Gt::Site_2::Weight                    Weight;

typedef CGAL::Apollonius_graph_2<Gt> Apollonius;

//typedef CGAL::Voronoi_diagram_2<Gt> VDA;

typedef Apollonius::Vertex_handle Vertex_handle;



using namespace std::chrono; 
int main(int argc, char **argv)
{

  /*std::ofstream myfile2;
  myfile2.open("FinalOutput2.txt");
	const int range_from  = 5000;
	const int range_to    = 10000;
	float num;
	for(int i=0;i<100000;i++)
	{
		std::random_device                  rand_dev;
		std::mt19937                        generator(rand_dev());
		std::uniform_int_distribution<int>  distr(range_from, range_to);
		num = distr(generator);
		myfile2<<num<<"\n";
	}
*/

	  std::ofstream myfile2;
  myfile2.open("FinalOutput2.txt");
	const int range_from  = -8500;
	const int range_to    = -5500;
	const int range_from1  = 2000;
	const int range_to1    = 6200;
	float num1,num2,ans;
	for(int i=0;i<1000;i++)
	{
		std::random_device                  rand_dev;
		std::mt19937                        generator(rand_dev());
		std::uniform_int_distribution<int>  distr(range_from, range_to);
		ans = distr(generator);
		num1 = ans/100;

		std::uniform_int_distribution<int>  distr1(range_from1, range_to1);
		ans = distr1(generator);
		num2 = ans/100;
		myfile2<<num1<<" "<<num2<<"\n";
	}

	return 0;
}