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

using namespace std;



/* Data Members of the structure Test
   x = x-coordinate of the point,
   y = y-coordinate of the point,
   w = weight of the point,
   euclid = Stores the Euclidean Distance between the point and Query Point
 */
struct Test 
{ 
   double x, y, w;
   //float euclid; 
};



/* A binary tree node has following data members
   data = Vector to store points,
   ag = stores additive vornoi diagram of the node,
   left = pointer to left child,
   right = pointer to right child */
struct Node { 
    //vector<Test> point;
    float px,py,pw;
    double maxprice,minprice;
    
};

struct Grid{
  double lx,ly,hx,hy;
  list<Node> g;
  
};

bool compare(Node a,Node b)
{
  if(a.maxprice > b.maxprice)
    return 1;
  else
    return 0;
}

using namespace std::chrono; 
int main(int argc, char **argv)
{
	//****PREPROCESSING  START*******
  //sleep(60);
  cout<<"PREPROCESSING STARTED"<<"\n";
  //sleep(30);
	//Open input file 
	std::ifstream ifs("sites.cin");
	assert( ifs );
	//Apollonius ag;      //Apol 
	Apollonius_site_2  site;
	std::vector<Test> val;
  while ( ifs >> site )
  {
    //ag.insert(site)

    val.push_back({site.x(),site.y(),site.weight()});
  }
  float min = 10000000;
  float max = -1000000;
  for(int i = 0;i<val.size();i++)
  {
      if(val[i].x>max)
        max = val[i].x;
      if(val[i].x<min)
        min = val[i].x;
  }

  cout<<max<<" "<<min<<"\n";

  min = 10000000;
  max = -1000000;
  for(int i = 0;i<val.size();i++)
  {
      if(val[i].y>max)
        max = val[i].y;
      if(val[i].y<min)
        min = val[i].y;
  }
   cout<<max<<" "<<min<<"\n";
}