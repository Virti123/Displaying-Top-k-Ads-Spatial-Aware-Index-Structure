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
#include <time.h>
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


bool compare(Test a,Test b)
{
  if(a.euclid<=b.euclid)
    return 1;
  else
    return 0;
}


  

int partition (std::vector<Test> &arr, int low, int high)  
{  
    float pivot = arr[high].euclid; // pivot  
    int i = (low - 1); // Index of smaller element  
    float temp1,temp2,temp3,temp4;
    for (int j = low; j <= high - 1; j++)  
    {  
        // If current element is smaller than the pivot  
        if (arr[j].euclid < pivot)  
        {  
            i++; // increment index of smaller element  
            //swap(&arr[i], &arr[j]);  
            temp1 = arr[i].x;
            temp2 = arr[i].y;
            temp3 = arr[i].w;
            temp4 = arr[i].euclid;
            arr[i].x = arr[j].x;
            arr[i].y = arr[j].y;
            arr[i].w = arr[j].w;
            arr[i].euclid = arr[j].euclid;
            arr[j].x = temp1;
            arr[j].y = temp2;
            arr[j].w = temp3;
            arr[j].euclid = temp4;
        }  
    }  
    //swap(&arr[i + 1], &arr[high]);  
    temp1 = arr[i+1].x;
    temp2 = arr[i+1].y;
    temp3 = arr[i+1].w;
    temp4 = arr[i+1].euclid;
    arr[i+1].x = arr[high].x;
    arr[i+1].y = arr[high].y;
    arr[i+1].w = arr[high].w;
    arr[i+1].euclid = arr[high].euclid;
    arr[high].x = temp1;
    arr[high].y = temp2;
    arr[high].w = temp3;
    arr[high].euclid = temp4;
    return (i + 1);  
}  
  

void quickSort(std::vector<Test> &arr, int low, int high)  
{  
    if (low < high)  
    {  
        /* pi is partitioning index, arr[p] is now  
        at right place */
        int pi = partition(arr, low, high);  
  
        // Separately sort elements before  
        // partition and after partition  
        quickSort(arr, low, pi - 1);  
        quickSort(arr, pi + 1, high);  
    }  
} 

using namespace std::chrono; 
int main(int argc, char **argv)
{
 // sleep(5);
  //Open input file 
  Apollonius_site_2  site;
  std::vector<Test> val;
  float a,b,euclidianDist;
  Point_2 p;
  //time_t start, end;
  int k=1;
  std::ifstream ips("queries1.dt.cin");
  std::ofstream myfile;
  std::ofstream timefile;
  myfile.open ("FinalOutput1.txt");
  timefile.open ("Time1.txt");
  while(ips>>p)
  {
    auto start_clock30 = high_resolution_clock::now(); 
    //clock_t tStart = clock();
    //auto start1 = high_resolution_clock::now(); 
    a = p.x();
    b = p.y(); 
  
  
  std::ifstream ifs("sites.cin");
  assert( ifs );
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
  quickSort(val, 0, val.size()-1);
  //for(int z=0;z<val.size();z++)
  //{
    //std::cout<<val[z].euclid<<"\n";
  //}

 myfile<<"-----------------------------------------------------"<<"\n";
  myfile<<"****The top "<<k<<" ads for points "<<a<<" "<<b<<" are **** " <<"\n";
  for(int z=0;z<k;z++)
  {
     myfile<<z+1<<")"<<val[z].x<<" "<<val[z].y<<"\n";
  }
  // auto stop = high_resolution_clock::now();
  // auto duration = duration_cast<microseconds>(stop - start1); 
  //std::cout << "\n\nTime taken by function: "<< duration.count() << " microseconds" ;
  //std::cout<<"\n";
  //printf("Time taken: %.7fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  auto stop_clock30 = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop_clock30 - start_clock30); 
  
  myfile << "\n\nTime taken by function: "
        << duration.count() << " microseconds"; 
  myfile<<"\n";
 timefile<< duration.count(); 
   timefile<<"\n";
 }
myfile.close();
timefile.close();
 ips.close();
  return 0;
}