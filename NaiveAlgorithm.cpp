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
#include<limits>
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


struct Test 
{ 
   double x, y, w,euclid; 
};


using namespace std::chrono; 
int main(int argc, char **argv)
{
 // sleep(5);
  //Open input file 
  Apollonius_site_2  site;
  std::vector<Test> val;
  float a,b,euclidianDist;
  Point_2 p;
   time_t start, end;
  std::ifstream ips("queries1.dt.cin");
  while(ips>>p)
  {
    a = p.x();
    b = p.y(); 
  
  
  std::ifstream ifs("sites.cin");
  assert( ifs );
  auto start1 = high_resolution_clock::now(); 
    time(&start);
    //ios_base::sync_with_stdio(false);
  //Apollonius ag;      //Apol 
  
  //create a vector from input file
  while ( ifs >> site )
  {
    euclidianDist = sqrt(pow( site.x()-a,2)+pow(site.y()-b,2))-site.weight();
    //ag.insert(site);
    val.push_back({site.x(),site.y(),site.weight(),euclidianDist});
  }
  ifs.close();
  int k=10;
  int j=0;

  while(k>0)
  {
    double min = std::numeric_limits<double>::infinity();
    int x=0;
    for(int i=0;i<val.size();i++)
    {
        if(val[i].euclid<min)
        {
          min = val[i].euclid;
          x=i;
        }
    }
    std::cout<<j<<")"<<val[x].x<<" "<<val[x].y<<"\n";
    j++;
    val.erase(val.begin() +x);
    k--;
}
time(&end); 
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start1); 
  
    std::cout << "\n\nTime taken by function: "
         << duration.count() << " microseconds"; 
    std::cout<<"\n";
    double time_taken = double(end - start); 
    std::cout << "Time taken by program is : "<< time_taken ; 
    std::cout << " sec "; 
    std::cout<<"\n\n";
 }
 ips.close();
  return 0;
}