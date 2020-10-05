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
   float x, y, w;
   //float euclid; 
};
struct OutputVector 
{ 
   float x, y, w, priceoffered;
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
    float maxprice,minprice;
    
};

struct Grid{
  float lx,ly,hx,hy;
  list<Node> g;
  
};

bool compare(Node a,Node b)
{
  if(a.maxprice > b.maxprice)
    return 1;
  else
    return 0;
}
bool compare1(OutputVector a,OutputVector b)
{
  if(a.priceoffered > b.priceoffered)
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
  Point_2 p1;

  while ( ifs >> site )
  {
    //ag.insert(site)

    val.push_back({site.x(),site.y(),site.weight()});
  }
  ifs.close();
  cout<<"Input read \n";
  //sleep(60);
  cout<<"start again \n";
  //cout<<sizeof(val[0]);
  int inputsize = val.size();
  int pcount=0;
  int gridsize=10;
  int gridsize_1 = gridsize-1;
  //cout<<"n = "<<n;
  Grid preprocess_data[gridsize][gridsize];
  Grid temp1[1][1];
  float xmin = -121.933;
  float xmax = 43.9118;

  float ymin = -75.5246;  
  float ymax = 51.8811;
  float p = xmin, q=ymin,m=0,n1=0;
  float a = (xmax - xmin)/gridsize;
  cout<<"a = "<<a<<"\n"; 
  float b = (ymax - ymin)/gridsize;
  
  float mindist=1000000,maxdist=0;
  float euclid = 0;
  for(int i=0;i<gridsize;i++)
  {
    p=xmin;
    m=0;
    q = q+n1;
    n1=b;
    for(int j=0;j<gridsize;j++)
    {
      preprocess_data[i][j].lx = p+m;
      p = preprocess_data[i][j].lx;
      preprocess_data[i][j].hx = preprocess_data[i][j].lx + a;
      m=a;
      preprocess_data[i][j].ly = q;
      preprocess_data[i][j].hy = preprocess_data[i][j].ly + b;
      //q = preprocess_data[i][j].ly;
      //cout<<preprocess_data[i][j].hx<<" "<<preprocess_data[i][j].hy<<" "<<" || ";
      //cout<<"value of n = "<<n;
      for(int k=0;k<inputsize;k++)
      {
        Node temp;
        mindist=1000000,maxdist=0;
        temp.px = val[k].x;
        temp.py = val[k].y;
        temp.pw = val[k].w;
        euclid = sqrt(pow(preprocess_data[i][j].lx-temp.px,2)+pow(preprocess_data[i][j].ly-temp.py,2));
        if(euclid>maxdist)
          maxdist=euclid;
        if(euclid<mindist)
          mindist=euclid;
        euclid = sqrt(pow(preprocess_data[i][j].lx-temp.px,2)+pow(preprocess_data[i][j].hy-temp.py,2));
         if(euclid>maxdist)
          maxdist=euclid;
        if(euclid<mindist)
          mindist=euclid;
        euclid = sqrt(pow(preprocess_data[i][j].hx-temp.px,2)+pow(preprocess_data[i][j].ly-temp.py,2));
         if(euclid>maxdist)
          maxdist=euclid;
        if(euclid<mindist)
          mindist=euclid;
        euclid = sqrt(pow(preprocess_data[i][j].hx-temp.px,2)+pow(preprocess_data[i][j].hy-temp.py,2));
         if(euclid>maxdist)
          maxdist=euclid;
        if(euclid<mindist)
          mindist=euclid;
        //cout<<"maxdist = "<<maxdist<<"\n";
        temp.maxprice = temp.pw - mindist;
        temp.minprice = temp.pw - maxdist;
        preprocess_data[i][j].g.push_back(temp);
          //temp1[0][0].g.push_back(temp);
        //cout<<preprocess_data[i][j].g[k].px<<" "<<preprocess_data[i][j].g[k].py<<" "<<" || ";
    
     }

 //preprocess_data[i][j].g.sort();
 preprocess_data[i][j].g.sort(compare);
    //   temp1[0][0].g.sort(compare);
// list<Node>::const_iterator first = temp1[0][0].g.begin();
//list<Node>::const_iterator last = temp1[0][0].g.begin() + 1000;
//preprocess_data[i][j].g.assign(std::next(temp1[0][0].g.begin(), 0), std::next(temp1[0][0].g.begin() , 1000));
// cout<<"Sizeeeee : "<<preprocess_data[i][j].g.size()<<"\n";


    }
  
    cout<<"\n";
    pcount++;
    cout<<"Hello"<<pcount<<"\n";
    //sleep(10);
    cout<<"start \n";
  }
  cout<<"Grid Created\n";
  //cout<<"Sizeeeee : "<<preprocess_data[0][0].g.size()<<"\n";
  //sleep(30);
  /*float count=0;
  list<Node>::iterator it;
  cout<<preprocess_data[0][0].lx<<" "<<preprocess_data[0][0].ly<<" "<<" || ";
   cout<<preprocess_data[0][0].hx<<" "<<preprocess_data[0][0].hy<<" "<<" || "<<"\n";
        it = preprocess_data[0][0].g.begin();
      while (it != preprocess_data[0][0].g.end()) {
        count++;
      cout << "maxprice: " << it->maxprice << " minprice: " << it->minprice <<"\n";
      it++;
      }
      cout<<"\n"<<count<<"\n";
return 0;*/
//cout<<"Hieeeeeeeeeeeeeeeeeeeee";
   int k=5;
  ofstream myfile2;
  ofstream timefile2;
  myfile2.open("FinalOutput2.txt");
  timefile2.open ("Time2.txt");
  std::ifstream ils("queries.dt.cin");
  assert( ils );
  float xval,yval;
  float distance=0;
  list<Node>::iterator it;
  list<Node>::iterator it1;
 
  vector<OutputVector> finaloutput;
  int flag=0;
  while(ils>>p1)
  {
    flag=0;
    auto start1 = high_resolution_clock::now(); 
    ios_base::sync_with_stdio(false);
    xval = p1.x();
    yval = p1.y();
    int countval=0;
    //cout<<"Hi Start\n";
    if(xval>=preprocess_data[0][0].lx && xval<=preprocess_data[gridsize_1][0].hx && yval>=preprocess_data[0][0].ly && yval<=preprocess_data[0][gridsize_1].hy)
    {
      cout<<"1\n";
      for(int i=0;i<gridsize_1;i++)
    {
      //cout<<preprocess_data[i][0].ly<<"<"<<yval<<"  "<<yval<<"<"<<preprocess_data[i][0].hy<<"\n";
      if(preprocess_data[i][0].ly<=yval && yval<=preprocess_data[i][0].hy)
      {
        //cout<<"Hi friend chail pe loo"<<"\n";
        for(int j=0;j<gridsize;j++)
        {
           //cout<<"fir se udd chala"<<"\n";
          if(preprocess_data[i][j].lx<xval && xval<preprocess_data[i][j].hx)
          {
            //cout<<"chal lonavala"<<"\n";

              it = preprocess_data[i][j].g.begin();
              it1 = it++;
               while (it1 != preprocess_data[i][j].g.end()) {
                //cout<<it->minprice<<">"<<it1->maxprice<<"  "<<countval<<">="<<k<<"\n";
                if(it->minprice>it1->maxprice && countval>=k)
                {
                  //cout<<"chalooooo lonavala 1111";
                  float maxvalue=-9999999;
                  int index=-1;
                  //sort(finaloutput.begin(),finaloutput.end(),compare1);
                  //cout<<"Size of Output = "<<finaloutput.size()<<"\n";
                  myfile2<<"-----------------------------------------------------------\n";
                  myfile2<<"*** Top "<<k<<" ads for points "<<xval<<" "<<yval<<" are ***"<<"\n";
                  for(int o=0;o<k;o++)
                  {
                    
                       index=-1;
                    for(int s=0;s<finaloutput.size();s++)
                    {
                      //cout<<"sssss";
                     
                        if (finaloutput[s].priceoffered>maxvalue)
                        {
                          maxvalue = finaloutput[s].priceoffered;
                          index=s;
                        }
                        
                    }
                    myfile2<<finaloutput[index].x<<" "<<finaloutput[index].y<<"\n";
                    finaloutput[index].priceoffered = -99999999;
                    maxvalue=-9999999;
                    
                  }
                  finaloutput.clear();
                  flag = 1;
                  break;
                }       
                else
                {
                  countval++;
                  distance = sqrt(pow(xval-it->px,2)+pow(yval-it->py,2));
                  float price = it->pw - distance;
                  finaloutput.push_back({it->px,it->py,it->pw,price});

                }
                it++;
                it1++;
          }
          }
          if(flag==1)
            break;
        }
      }
      if(flag==1)
        break;
    }   
    }
    else if(yval<preprocess_data[0][0].ly && xval>=preprocess_data[0][0].lx && xval<=preprocess_data[gridsize_1][0].hx)
    {
      //check all data[0][0] to data[29][0] x values
      cout<<"2\n";
      float x1,y1,d;
      int r;
      float min = 9999999;
      y1 = preprocess_data[0][0].ly;
       for(int j=0;j<gridsize;j++)
        {
           //cout<<"fir se udd chala"<<"\n";
          d = sqrt(pow(xval-preprocess_data[j][0].lx,2)+pow(yval-y1,2));
          if(d<min)
          {
            min=d;
            r=j;
          }
          
        }
          

            //cout<<"chal lonavala"<<"\n";

              it = preprocess_data[r][0].g.begin();
              it1 = it++;
               while (it1 != preprocess_data[r][0].g.end()) {
                //cout<<it->minprice<<">"<<it1->maxprice<<"  "<<countval<<">="<<k<<"\n";
                if(it->minprice>it1->maxprice && countval>=k)
                {
                  //cout<<"chalooooo lonavala 1111";
                  float maxvalue=-9999999;
                  int index=-1;
                  //sort(finaloutput.begin(),finaloutput.end(),compare1);
                  //cout<<"Size of Output = "<<finaloutput.size()<<"\n";
                   myfile2<<"-----------------------------------------------------------\n";
                  myfile2<<"*** Top "<<k<<" ads for points "<<xval<<" "<<yval<<" are ***"<<"\n";
                  for(int o=0;o<k;o++)
                  {
                    
                       index=-1;
                    for(int s=0;s<finaloutput.size();s++)
                    {
                      //cout<<"sssss";
                     
                        if (finaloutput[s].priceoffered>maxvalue)
                        {
                          maxvalue = finaloutput[s].priceoffered;
                          index=s;
                        }
                        
                    }
                        myfile2<<finaloutput[index].x<<" "<<finaloutput[index].y<<" "<< "\n";
                        finaloutput[index].priceoffered = -99999999;
                        maxvalue=-9999999;
                    
                  }
                  finaloutput.clear();
                  flag = 1;
                  break;
                }       
                else
                {
                  countval++;
                  distance = sqrt(pow(xval-it->px,2)+pow(yval-it->py,2));
                  float price = it->pw - distance;
                  finaloutput.push_back({it->px,it->py,it->pw,price});

                }
                it++;
                it1++;
          }
        
    }
    else if(yval>=preprocess_data[0][gridsize_1].hy && xval>=preprocess_data[0][0].lx && xval<=preprocess_data[gridsize_1][0].hx)
    {
      //check all data[0][29] to data[29][29] x values
      cout<<"3\n";
      float x1,y1,d;
      int r;
      float min = 9999999;
      y1 = preprocess_data[0][gridsize_1].ly;
       for(int j=0;j<gridsize;j++)
        {
           //cout<<"fir se udd chala"<<"\n";
          d = sqrt(pow(xval-preprocess_data[j][gridsize_1].lx,2)+pow(yval-y1,2));
          if(d<min)
          {
            min=d;
            r=j;
          }
          
        }
          

            //cout<<"chal lonavala"<<"\n";

              it = preprocess_data[r][gridsize_1].g.begin();
              it1 = it++;
               while (it1 != preprocess_data[r][gridsize_1].g.end()) {
                //cout<<it->minprice<<">"<<it1->maxprice<<"  "<<countval<<">="<<k<<"\n";
                if(it->minprice>it1->maxprice && countval>=k)
                {
                  //cout<<"chalooooo lonavala 1111";
                  float maxvalue=-9999999;
                  int index=-1;
                  //sort(finaloutput.begin(),finaloutput.end(),compare1);
                  //cout<<"Size of Output = "<<finaloutput.size()<<"\n";
                   myfile2<<"-----------------------------------------------------------\n";
                  myfile2<<"*** Top "<<k<<" ads for points "<<xval<<" "<<yval<<" are ***"<<"\n";
                  for(int o=0;o<k;o++)
                  {
                    
                       index=-1;
                    for(int s=0;s<finaloutput.size();s++)
                    {
                      //cout<<"sssss";
                     
                        if (finaloutput[s].priceoffered>maxvalue)
                        {
                          maxvalue = finaloutput[s].priceoffered;
                          index=s;
                        }
                        
                    }
                        myfile2<<finaloutput[index].x<<" "<<finaloutput[index].y<<"\n";
                        finaloutput[index].priceoffered = -99999999;
                        maxvalue=-9999999;
                    
                  }
                  finaloutput.clear();
                  flag = 1;
                  break;
                }       
                else
                {
                  countval++;
                  distance = sqrt(pow(xval-it->px,2)+pow(yval-it->py,2));
                  float price = it->pw - distance;
                  finaloutput.push_back({it->px,it->py,it->pw,price});

                }
                it++;
                it1++;
          }
        
    }
    else if(xval<=preprocess_data[0][0].lx && yval>=preprocess_data[0][0].ly && yval<=preprocess_data[0][gridsize_1].hy)
    {
      cout<<"4\n";
      //check all data[0][0] to data[0][29] y values
      float x1,y1,d;
      int r;
      float min = 9999999;
      x1 = preprocess_data[0][0].lx;
       for(int j=0;j<gridsize;j++)
        {
           //cout<<"fir se udd chala"<<"\n";
          d = sqrt(pow(xval-x1,2)+pow(yval-preprocess_data[0][j].ly,2));
          if(d<min)
          {
            min=d;
            r=j;
          }
          
        }
          

            //cout<<"chal lonavala"<<"\n";

              it = preprocess_data[0][r].g.begin();
              it1 = it++;
               while (it1 != preprocess_data[0][r].g.end()) {
                //cout<<it->minprice<<">"<<it1->maxprice<<"  "<<countval<<">="<<k<<"\n";
                if(it->minprice>it1->maxprice && countval>=k)
                {
                  //cout<<"chalooooo lonavala 1111";
                  float maxvalue=-9999999;
                  int index=-1;
                  //sort(finaloutput.begin(),finaloutput.end(),compare1);
                  //cout<<"Size of Output = "<<finaloutput.size()<<"\n";
                   myfile2<<"-----------------------------------------------------------\n";
                  myfile2<<"*** Top "<<k<<" ads for points "<<xval<<" "<<yval<<" are ***"<<"\n";
                  for(int o=0;o<k;o++)
                  {
                    
                       index=-1;
                    for(int s=0;s<finaloutput.size();s++)
                    {
                      //cout<<"sssss";
                     
                        if (finaloutput[s].priceoffered>maxvalue)
                        {
                          maxvalue = finaloutput[s].priceoffered;
                          index=s;
                        }
                        
                    }
                        myfile2<<finaloutput[index].x<<" "<<finaloutput[index].y<<"\n";
                        finaloutput[index].priceoffered = -99999999;
                        maxvalue=-9999999;
                    
                  }
                  finaloutput.clear();
                  flag = 1;
                  break;
                }       
                else
                {
                  countval++;
                  distance = sqrt(pow(xval-it->px,2)+pow(yval-it->py,2));
                  float price = it->pw - distance;
                  finaloutput.push_back({it->px,it->py,it->pw,price});

                }
                it++;
                it1++;
          }
    }
    else if(xval>=preprocess_data[0][gridsize_1].lx && yval>=preprocess_data[0][0].ly && yval<=preprocess_data[gridsize_1][gridsize_1].hy)
    {
        //check all data[29][0] to data[29][29] y values
      cout<<"5\n";
      float x1,y1,d;
      int r;
      float min = 9999999;
      x1 = preprocess_data[gridsize_1][0].lx;
       for(int j=0;j<gridsize;j++)
        {
           //cout<<"fir se udd chala"<<"\n";
          d = sqrt(pow(xval-x1,2)+pow(yval-preprocess_data[gridsize_1][j].ly,2));
          if(d<min)
          {
            min=d;
            r=j;
          }
          
        }
          

            //cout<<"chal lonavala"<<"\n";

              it = preprocess_data[gridsize_1][r].g.begin();
              it1 = it++;
               while (it1 != preprocess_data[gridsize_1][r].g.end()) {
                //cout<<it->minprice<<">"<<it1->maxprice<<"  "<<countval<<">="<<k<<"\n";
                if(it->minprice>it1->maxprice && countval>=k)
                {
                  //cout<<"chalooooo lonavala 1111";
                  float maxvalue=-9999999;
                  int index=-1;
                  //sort(finaloutput.begin(),finaloutput.end(),compare1);
                  //cout<<"Size of Output = "<<finaloutput.size()<<"\n";
                   myfile2<<"-----------------------------------------------------------\n";
                  myfile2<<"*** Top "<<k<<" ads for points "<<xval<<" "<<yval<<" are ***"<<"\n";
                  for(int o=0;o<k;o++)
                  {
                    
                       index=-1;
                    for(int s=0;s<finaloutput.size();s++)
                    {
                      //cout<<"sssss";
                     
                        if (finaloutput[s].priceoffered>maxvalue)
                        {
                          maxvalue = finaloutput[s].priceoffered;
                          index=s;
                        }
                        
                    }
                        myfile2<<finaloutput[index].x<<" "<<finaloutput[index].y<<"\n";
                        finaloutput[index].priceoffered = -99999999;
                        maxvalue=-9999999;
                    
                  }
                  finaloutput.clear();
                  flag = 1;
                  break;
                }       
                else
                {
                  countval++;
                  distance = sqrt(pow(xval-it->px,2)+pow(yval-it->py,2));
                  float price = it->pw - distance;
                  finaloutput.push_back({it->px,it->py,it->pw,price});

                }
                it++;
                it1++;
          }

    }    
    else if(yval<preprocess_data[0][0].ly && xval<preprocess_data[0][0].lx)
    {

      //only check preprocess_data[0][0]
      cout<<"6\n";

              it = preprocess_data[0][0].g.begin();
              it1 = it++;
               while (it1 != preprocess_data[0][0].g.end()) {
                //cout<<it->minprice<<">"<<it1->maxprice<<"  "<<countval<<">="<<k<<"\n";
                if(it->minprice>it1->maxprice && countval>=k)
                {
                  //cout<<"chalooooo lonavala 1111";
                  float maxvalue=-9999999;
                  int index=-1;
                  //sort(finaloutput.begin(),finaloutput.end(),compare1);
                  //cout<<"Size of Output = "<<finaloutput.size()<<"\n";
                   myfile2<<"-----------------------------------------------------------\n";
                  myfile2<<"*** Top "<<k<<" ads for points "<<xval<<" "<<yval<<" are ***"<<"\n";
                  for(int o=0;o<k;o++)
                  {
                    
                       index=-1;
                    for(int s=0;s<finaloutput.size();s++)
                    {
                      //cout<<"sssss";
                     
                        if (finaloutput[s].priceoffered>maxvalue)
                        {
                          maxvalue = finaloutput[s].priceoffered;
                          index=s;
                        }
                        
                    }
                    myfile2<<finaloutput[index].x<<" "<<finaloutput[index].y<<"\n";
                    finaloutput[index].priceoffered = -99999999;
                    maxvalue=-9999999;
                    
                  }
                  finaloutput.clear();
                  flag = 1;
                  break;
                }       
                else
                {
                  countval++;
                  distance = sqrt(pow(xval-it->px,2)+pow(yval-it->py,2));
                  float price = it->pw - distance;
                  finaloutput.push_back({it->px,it->py,it->pw,price});

                }
                it++;
                it1++;
          }


    }
    else if(xval>=preprocess_data[gridsize_1][0].hx && yval<=preprocess_data[0][0].ly)
    {
      //check only data[29][0]
      cout<<"7\n";

              it = preprocess_data[gridsize_1][0].g.begin();
              it1 = it++;
               while (it1 != preprocess_data[gridsize_1][0].g.end()) {
                //cout<<it->minprice<<">"<<it1->maxprice<<"  "<<countval<<">="<<k<<"\n";
                if(it->minprice>it1->maxprice && countval>=k)
                {
                  //cout<<"chalooooo lonavala 1111";
                  float maxvalue=-9999999;
                  int index=-1;
                  //sort(finaloutput.begin(),finaloutput.end(),compare1);
                  //cout<<"Size of Output = "<<finaloutput.size()<<"\n";
                   myfile2<<"-----------------------------------------------------------\n";
                  myfile2<<"*** Top "<<k<<" ads for points "<<xval<<" "<<yval<<" are ***"<<"\n";
                  for(int o=0;o<k;o++)
                  {
                    
                       index=-1;
                    for(int s=0;s<finaloutput.size();s++)
                    {
                      //cout<<"sssss";
                     
                        if (finaloutput[s].priceoffered>maxvalue)
                        {
                          maxvalue = finaloutput[s].priceoffered;
                          index=s;
                        }
                        
                    }
                    myfile2<<finaloutput[index].x<<" "<<finaloutput[index].y<<"\n";
                    finaloutput[index].priceoffered = -99999999;
                    maxvalue=-9999999;
                    
                  }
                  finaloutput.clear();
                  flag = 1;
                  break;
                }       
                else
                {
                  countval++;
                  distance = sqrt(pow(xval-it->px,2)+pow(yval-it->py,2));
                  float price = it->pw - distance;
                  finaloutput.push_back({it->px,it->py,it->pw,price});

                }
                it++;
                it1++;
          }
    }
    else if(xval<=preprocess_data[0][gridsize_1].lx && yval>=preprocess_data[0][gridsize_1].hy)
    {
        //check only data[0][29]
          cout<<"8\n";

              it = preprocess_data[0][gridsize_1].g.begin();
              it1 = it++;
               while (it1 != preprocess_data[0][gridsize_1].g.end()) {
                //cout<<"p\n";
                //cout<<it->minprice<<">"<<it1->maxprice<<"  "<<countval<<">="<<k<<"\n";
                if(it->minprice>it1->maxprice && countval>=k)
                {
                  
                  break;
                }       
                else
                {
                 // cout<<"In else\n";
                  countval++;
                  distance = sqrt(pow(xval-it->px,2)+pow(yval-it->py,2));
                  float price = it->pw - distance;
                  finaloutput.push_back({it->px,it->py,it->pw,price});

                }
                it++;
                it1++;
          }
         // cout<<"In it<minprice and countval>=k\n";
                  //cout<<"chalooooo lonavala 1111";
                  float maxvalue=-9999999;
                  int index=-1;
                  //sort(finaloutput.begin(),finaloutput.end(),compare1);
                 // cout<<"Size of Output = "<<finaloutput.size()<<"\n";
                   myfile2<<"-----------------------------------------------------------\n";
                  myfile2<<"*** Top "<<k<<" ads for points "<<xval<<" "<<yval<<" are ***"<<"\n";
                  for(int o=0;o<k;o++)
                  {
                    
                       index=-1;
                    for(int s=0;s<finaloutput.size();s++)
                    {
                      //cout<<"sssss";
                     
                        if (finaloutput[s].priceoffered>maxvalue)
                        {
                          maxvalue = finaloutput[s].priceoffered;
                          index=s;
                        }
                        
                    }
                    myfile2<<finaloutput[index].x<<" "<<finaloutput[index].y<<"\n";
                    finaloutput[index].priceoffered = -99999999;
                    maxvalue=-9999999;
                    
                  }
                  finaloutput.clear();
                  flag = 1;
    }
    else if(xval>=preprocess_data[0][gridsize_1].hx && yval>=preprocess_data[gridsize_1][gridsize_1].hy)
    {
        //check only data[29][29]
      cout<<"9\n";

              it = preprocess_data[gridsize_1][gridsize_1].g.begin();
              it1 = it++;
               while (it1 != preprocess_data[gridsize_1][gridsize_1].g.end()) {
                //cout<<it->minprice<<">"<<it1->maxprice<<"  "<<countval<<">="<<k<<"\n";
                if(it->minprice>it1->maxprice && countval>=k)
                {
                  //cout<<"chalooooo lonavala 1111";
                  
                  break;
                }       
                else
                {
                  countval++;
                  distance = sqrt(pow(xval-it->px,2)+pow(yval-it->py,2));
                  float price = it->pw - distance;
                  finaloutput.push_back({it->px,it->py,it->pw,price});

                }
                it++;
                it1++;
          }
          float maxvalue=-9999999;
                  int index=-1;
                  //sort(finaloutput.begin(),finaloutput.end(),compare1);
                 // cout<<"Size of Output = "<<finaloutput.size()<<"\n";
                   myfile2<<"-----------------------------------------------------------\n";
                  myfile2<<"*** Top "<<k<<" ads for points "<<xval<<" "<<yval<<" are ***"<<"\n";
                  for(int o=0;o<k;o++)
                  {
                    
                       index=-1;
                    for(int s=0;s<finaloutput.size();s++)
                    {
                      //cout<<"sssss";
                     
                        if (finaloutput[s].priceoffered>maxvalue)
                        {
                          maxvalue = finaloutput[s].priceoffered;
                          index=s;
                        }
                        
                    }
                    myfile2<<finaloutput[index].x<<" "<<finaloutput[index].y<<"\n";
                    finaloutput[index].priceoffered = -99999999;
                    maxvalue=-9999999;
                    
                  }
                  finaloutput.clear();
                  flag = 1;
    }


    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start1); 
    myfile2<< "\n\nTime taken by function: "<< duration.count() << " microseconds" ;
    myfile2<<"\n";
    timefile2<<duration.count();
    timefile2<<"\n";

    
  }
/*std::ifstream iqs("queries.dt.cin");
  assert( iqs );
  //Apollonius ag;      //Apol 

  while ( iqs >> p1 )
  {
    //cout<<"Guess I am inn"<<"\n";
    //ag.insert(site)
    xval = p1.x();
    yval = p1.y();
    int countval=0;
    for(int i=0;i<29;i++)
    {
      //cout<<preprocess_data[i][0].ly<<"<"<<yval<<"  "<<yval<<"<"<<preprocess_data[i][0].hy<<"\n";
      if(preprocess_data[i][0].ly<=yval && yval<=preprocess_data[i][0].hy)
      {
        //cout<<"Hi friend chail pe loo"<<"\n";
        for(int j=0;j<30;j++)
        {
           //cout<<"fir se udd chala"<<"\n";
          if(preprocess_data[i][j].lx<xval && xval<preprocess_data[i][j].hx)
          {
            //cout<<"chal lonavala"<<"\n";

              it = preprocess_data[i][j].g.begin();
              it1 = it++;
               while (it1 != preprocess_data[i][j].g.end()) {
                //cout<<it->minprice<<">"<<it1->maxprice<<"  "<<countval<<">="<<k<<"\n";
                if(it->minprice>it1->maxprice && countval>=k)
                {
                  cout<<"chalooooo lonavala 1111";
                  float maxvalue=-9999999;
                  int index=-1;
                  //sort(finaloutput.begin(),finaloutput.end(),compare1);
                  cout<<"Size of Output = "<<finaloutput.size()<<"\n";
                  cout<<"*** Top "<<k<<" ads for points "<<xval<<" "<<yval<<" are ***"<<"\n";
                  for(int o=0;o<k;o++)
                  {
                    
                       index=-1;
                    for(int s=0;s<finaloutput.size();s++)
                    {
                      //cout<<"sssss";
                     
                        if (finaloutput[s].priceoffered>maxvalue)
                        {
                          maxvalue = finaloutput[s].priceoffered;
                          index=s;
                        }
                        
                    }
                    cout<<finaloutput[index].x<<" "<<finaloutput[index].y<<" "<< finaloutput[index].priceoffered<<"\n";
                        finaloutput[index].priceoffered = -99999999;
                        maxvalue=-9999999;
                    
                  }
                  finaloutput.clear();
                  flag = 1;
                  break;
                }       
                else
                {
                  countval++;
                  distance = sqrt(pow(xval-it->px,2)+pow(yval-it->py,2));
                  float price = it->pw - distance;
                  finaloutput.push_back({it->px,it->py,it->pw,price});

                }
                it++;
                it1++;
          }
          }
          if(flag==1)
            break;
        }
      }
      if(flag==1)
        break;
    }
    
  }*/


	return 0;
}