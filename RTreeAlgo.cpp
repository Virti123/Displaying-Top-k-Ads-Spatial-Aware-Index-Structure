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
#include<queue>
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
   // double euclid; 
};





/* A binary tree node has following data members
   data = Vector to store points,
   ag = stores additive vornoi diagram of the node,
   left = pointer to left child,
   right = pointer to right child */
struct Node { 
  vector<Test> point;
   double xmin, ymin, xmax,ymax,minweight,count,maxweight,minprice,maxprice;
  Node* pointer1;
  Node* pointer2;
  Node* pointer3;
  Node* pointer4;
    
};


struct compare { 
    bool operator()(Node* const& p1, Node* const& p2) 
    { 
        // return "true" if "p1" is ordered  
        // before "p2", for example: 
        return p1->minprice<p2->minprice; 
    } 
}; 

struct compare1 { 
    bool operator()(Test const& p1, Test const& p2) 
    { 
        // return "true" if "p1" is ordered  
        // before "p2", for example: 
        return p1.w<p2.w; 
    } 
}; 



void printPostorder(struct Node* node) 
{ 
    Point_2 p;
    Apollonius_site_2 site;
    if (node == NULL) 
        return; 
  
    // first recur on left subtree 
    if(node->pointer1==NULL && node->pointer2==NULL && node->pointer3==NULL && node->pointer4==NULL)
    {  
      
      //printf("(%lf,%lf,%lf) || ",p[i].x,p[i].y,p[i].w);
      cout<<"leaf\n";
     for(int a=0;a<node->point.size();a++)
     {
      cout<<node->point[a].x<<" "<<node->point[a].y<<" "<<node->point[a].w<<"\n";
      }
    
      //cout<<node->point[1].xmin<<" "<<node->point[1].ymin<<" "<<node->point[1].minweight<<"\n";
      //cout<<node->point[2].xmin<<" "<<node->point[2].ymin<<" "<<node->point[2].minweight<<"\n";
      //cout<<node->point[3].xmin<<" "<<node->point[3].ymin<<" "<<node->point[3].minweight<<"\n";
    }

    else{
       
      cout<<node->xmin<<" "<<node->ymin<<" "<<node->xmax<<" "<<node->ymax<<" "<<node->minweight<<" "<<node->maxweight<<" "<<node->count<<"\n";
    
    //cout<<node->point[0].xmin<<" "<<node->point[0].ymin<<" "<<node->point[0].xmax<<" "<<node->point[0].ymax<<"\n";
    //cout<<node->point[1].xmin<<" "<<node->point[1].ymin<<" "<<node->point[1].xmax<<" "<<node->point[0].ymax<<"\n";
    //cout<<node->point[2].xmin<<" "<<node->point[2].ymin<<" "<<node->point[2].xmax<<" "<<node->point[0].ymax<<"\n";
    //cout<<node->point[3].xmin<<" "<<node->point[3].ymin<<" "<<node->point[3].xmax<<" "<<node->point[0].ymax<<"\n";
      cout<<"Pointer1\n";
    printPostorder(node->pointer1);
    cout<<"Pointer2\n";
    printPostorder(node->pointer2);
    cout<<"Pointer3\n";
    printPostorder(node->pointer3); 
    cout<<"Pointer4\n";
    printPostorder(node->pointer4);
  }
  
    // then recur on right subtree
    // now deal with the node. if the node is leaf node
    //print all the members of data vector of tree node
   
      printf("\n\n");
}



vector<Test> temp;
struct Node* createRTree(vector<Test> &v,double x1,double y1,double x2,double y2)
{
  
  //cout<<"---------------"<<"\n";
  int count=0;
  int c=0;
  struct Node* newnode;
  for(int i=0;i<v.size();i++)
  {
    if(x1<=v[i].x && v[i].x<=x2 && y1<=v[i].y && v[i].y<=y2)
    {
      count++;
      temp.push_back({v[i].x,v[i].y,v[i].w});
      if(count>4)
      {
        break;
      }
    }
  }
  if(count==0)
  {
   // cout<<"In NULL";
    return NULL;
  }
      if(count>4)
      {
        //cout<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<"\n";
        temp.clear();
       // list<Node> listNode;
       newnode = new Node(); 
       
        // double p1 = x1, q1 = (y1+y2)/2, r1 = (x1+x2)/2, s1 = y2,
        double max = -9999999,min=99999999;
        for(int j=0;j<v.size();j++)
        {
          if(x1<=v[j].x && v[j].x<=x2 && y1<=v[j].y && v[j].y<=y2)
          {
            c++;
            if(v[j].w>max)
            {
              max = v[j].w;
            }
            if(v[j].w<min)
            {
              min=v[j].w;
            }
          }
        }
        //cout<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<" "<<c<<"\n";
        newnode->xmin = x1;newnode->ymin = y1;newnode->xmax = x2;newnode->ymax = y2;newnode->minweight=min;newnode->maxweight=max;
        newnode->count=c;newnode->minprice=0;newnode->maxprice=0;
        newnode->pointer1 = createRTree(v,x1,(y1+y2)/2,(x1+x2)/2,y2);
        newnode->pointer2 = createRTree(v,(x1+x2)/2,(y1+y2)/2,x2,y2);
        newnode->pointer3 = createRTree(v,x1,y1,(x1+x2)/2,(y1+y2)/2);
        newnode->pointer4 = createRTree(v,(x1+x2)/2,y1,x2,(y1+y2)/2);
      
    }
  
  if(count<=4)
  {
     //cout<<"In count less than 4";
     //cout<<"Leaf Node\n";
      newnode = new Node();
      for(int a=0;a<temp.size();a++)
      {
        //cout<<temp[a].x<<" "<<temp[a].y<<" "<<temp[a].w<<"\n";

        newnode->point.push_back({temp[a].x,temp[a].y,temp[a].w});
      }
      newnode->count=-1;
      newnode->pointer1 = NULL;
      newnode->pointer2 = NULL;
      newnode->pointer3 = NULL;
      newnode->pointer4 = NULL;
      temp.clear();
  }


  return newnode;
}


void calculateMinMaxPrice(double xval,double yval,struct Node* node1)
{
  double euclid, mindist=1000000,maxdist=0;
  euclid = sqrt(pow(node1->xmin-xval,2)+pow(node1->ymin-yval,2));
        if(euclid>maxdist)
          maxdist=euclid;
        if(euclid<mindist)
          mindist=euclid;
        euclid = sqrt(pow(node1->xmax-xval,2)+pow(node1->ymin-yval,2));
         if(euclid>maxdist)
          maxdist=euclid;
        if(euclid<mindist)
          mindist=euclid;
        euclid = sqrt(pow(node1->xmin-xval,2)+pow(node1->ymax-yval,2));
         if(euclid>maxdist)
          maxdist=euclid;
        if(euclid<mindist)
          mindist=euclid;
        euclid = sqrt(pow(node1->xmax-xval,2)+pow(node1->ymax-yval,2));
         if(euclid>maxdist)
          maxdist=euclid;
        if(euclid<mindist)
          mindist=euclid;
        node1->minprice = node1->minweight - maxdist;
        node1->maxprice = node1->maxweight - mindist;

}

priority_queue<Test,vector<Test>,compare1> finalOutput;
double getOutput( priority_queue<Node*,vector<Node*>,compare> &q,double maxmin,double xcor,double ycor)
{
    struct Node* a = q.top();
    q.pop();
    // cout<<"Val = "<<a->minprice<<" "<<a->maxprice<<" "<<a->count<<"\n";
    //cout<<"Val = "<<a->minprice<<" "<<a->maxprice<<" "<<a->count<<"\n";

    if(a->count!=0)
    {
   
    priority_queue<Node*,vector<Node*>,compare> tempq;
    if(a->pointer1 == NULL && a->pointer2==NULL && a->pointer3==NULL && a->pointer4==NULL)
    {
      //cout<<"Printtttttttttttttting   MY NAME IS JOKER\n";
      for(int u=0;u<a->point.size();u++)
      {
      double euclid = sqrt(pow(a->point[u].x-xcor,2)+pow(a->point[u].y-ycor,2));
      double price = a->point[u].w - euclid;
      Test ans = {a->point[u].x,a->point[u].y,price};
      //cout<<a->point[u].x<<" "<<a->point[u].y<<" "<<price<<"\n";
      finalOutput.push(ans);
    }
      
    }
    else
    {
     // cout<<"\nIn else\n";
              //cout<<"\nIn if maxmin\n";
      if(a->maxprice>maxmin){
        if(a->pointer1!=NULL)
        {
          if((a->pointer1)->pointer1 == NULL && (a->pointer1)->pointer2==NULL && (a->pointer1)->pointer3==NULL && (a->pointer1)->pointer4==NULL)
          {
             // cout<<"Printtttttttttttttting MY NAME IS Vertex_handle\n";
              for(int u=0;u<(a->pointer1)->point.size();u++)
              {
              double euclid = sqrt(pow((a->pointer1)->point[u].x-xcor,2)+pow((a->pointer1)->point[u].y-ycor,2));
              double price = (a->pointer1)->point[u].w - euclid;
              Test ans = {(a->pointer1)->point[u].x,(a->pointer1)->point[u].y,price};
             // cout<<(a->pointer1)->point[u].x<<" "<<(a->pointer1)->point[u].y<<" "<<price<<"\n";
              finalOutput.push(ans);
              }
          }
          else
          {

          calculateMinMaxPrice(xcor,ycor,a->pointer1);
          if ((a->pointer1)->maxprice>maxmin)
            {
             tempq.push(a->pointer1);
            }
           
          }
        }

        if(a->pointer2!=NULL)
        {
          if((a->pointer2)->pointer1 == NULL && (a->pointer2)->pointer2==NULL && (a->pointer2)->pointer3==NULL && (a->pointer2)->pointer4==NULL)
          {
              //cout<<"Printtttttttttttttting\n";
              for(int u=0;u<(a->pointer2)->point.size();u++)
              {
              double euclid = sqrt(pow((a->pointer2)->point[u].x-xcor,2)+pow((a->pointer2)->point[u].y-ycor,2));
              double price = (a->pointer2)->point[u].w - euclid;
              Test ans = {(a->pointer2)->point[u].x,(a->pointer2)->point[u].y,price};
              //cout<<(a->pointer2)->point[u].x<<" "<<(a->pointer2)->point[u].y<<" "<<price<<"\n";
              finalOutput.push(ans);
            }
          }
          else
          {

          calculateMinMaxPrice(xcor,ycor,a->pointer2);
          if ((a->pointer2)->maxprice>maxmin)
            {
             tempq.push(a->pointer2);
            }
           
          }
        }
        

        if(a->pointer3!=NULL)
        {
          if((a->pointer3)->pointer1 == NULL && (a->pointer3)->pointer2==NULL && (a->pointer3)->pointer3==NULL && (a->pointer3)->pointer4==NULL)
          {
              //cout<<"Printtttttttttttttting\n";
              for(int u=0;u<(a->pointer3)->point.size();u++)
              {
              double euclid = sqrt(pow((a->pointer3)->point[u].x-xcor,2)+pow((a->pointer3)->point[u].y-ycor,2));
              double price = (a->pointer3)->point[u].w - euclid;
              Test ans = {(a->pointer3)->point[u].x,(a->pointer3)->point[u].y,price};
              //cout<<(a->pointer3)->point[u].x<<" "<<(a->pointer3)->point[u].y<<" "<<price<<"\n";
              finalOutput.push(ans);
            }
          }
          else
          {

          calculateMinMaxPrice(xcor,ycor,a->pointer3);
          if ((a->pointer3)->maxprice>maxmin)
            {
             tempq.push(a->pointer3);
            }
           
          }
        }

       if(a->pointer4!=NULL)
        {
          if((a->pointer4)->pointer1 == NULL && (a->pointer4)->pointer2==NULL && (a->pointer4)->pointer3==NULL && (a->pointer4)->pointer4==NULL)
          {
              //cout<<"Printtttttttttttttting\n";
              for(int u=0;u<(a->pointer4)->point.size();u++)
              {
              double euclid = sqrt(pow((a->pointer4)->point[u].x-xcor,2)+pow((a->pointer4)->point[u].y-ycor,2));
              double price = (a->pointer4)->point[u].w - euclid;
              Test ans = {(a->pointer4)->point[u].x,(a->pointer4)->point[u].y,price};
             // cout<<(a->pointer4)->point[u].x<<" "<<(a->pointer4)->point[u].y<<" "<<price<<"\n";
              finalOutput.push(ans);
            }
          }
          else
          {

          calculateMinMaxPrice(xcor,ycor,a->pointer4);
          if ((a->pointer4)->maxprice>maxmin)
            {
             tempq.push(a->pointer4);
            }
           
          }
        }
               
       
        while(!tempq.empty())
        {
          struct Node* s = tempq.top();
          
           // cout<<"Pushed";
         //   maxmin = s->minprice;
            q.push(s);
          
          tempq.pop();
        }
        
        //int countVal=0;
    }
  }
    
  }
  
  //cout<<"Queue Size = "<<q.size()<<"\n";
  struct Node* l = q.top();
  maxmin = l->minprice;
    return maxmin;

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
  double xmin1 = -121.933;
  double xmax1 = 43.9118;
  double ymin1 = -75.5246;  
  double ymax1 =51.8811 ;
  struct Node *head  = createRTree(val,xmin1,ymin1,xmax1,ymax1);
  cout<<"Doooooooooooooooooonnnnnnnnnnnnnnnneeeeeeeeeeee\n";
  priority_queue<Node*,vector<Node*>,compare> finalQueue;
 // printPostorder(head);


 std::ifstream ils("queries2.dt.cin");
  std::ofstream myfile2;
  std::ofstream timefile2;
  myfile2.open("FinalOutput2.txt");
  timefile2.open("TimeFile2.txt");
 // assert( ils );
  double xval,yval;
// cout<<"Helllllllllllllloooooooooooooo11111\n";
  while(ils>>p1)
  {
    int k=2;
     //cout<<"Helllllllllllllloooooooooooooo22222\n";
    auto start1 = high_resolution_clock::now(); 
    ios_base::sync_with_stdio(false);
    //cout<<"Helllllllllllllloooooooooooooo33333\n";
    xval=p1.x();
   // cout<<"Helllllllllllllloooooooooooooo44444\n";
    yval=p1.y(); 
   // cout<<"Helllllllllllllloooooooooooooo55555\n";
    myfile2<<"-----------------------------------------------------------------------\n";
    myfile2<<"*****The top-"<<k<<" ads for points"<<xval<<" "<<yval<<" are ****"<<"\n";
   // cout<<"-----------------------------------------------------------------------\n";
   // cout<<"*****The top-"<<k<<" ads for points"<<xval<<" "<<yval<<" are ****"<<"\n";
    priority_queue<Node*,vector<Node*>,compare> q;

    struct Node* treeNode = head;

    //struct Node* node1 = treeNode->pointer1;
    if(treeNode->pointer1!=NULL)
    {
    if((treeNode->pointer1)->pointer1==NULL && (treeNode->pointer1)->pointer2==NULL && (treeNode->pointer1)->pointer3==NULL && (treeNode->pointer1)->pointer4==NULL)
    {
         for(int u=0;u<(treeNode->pointer1)->point.size();u++)
      {
      double euclid = sqrt(pow((treeNode->pointer1)->point[u].x-xval,2)+pow((treeNode->pointer1)->point[u].y-yval,2));
      double price = (treeNode->pointer1)->point[u].w - euclid;
      Test ans = {(treeNode->pointer1)->point[u].x,(treeNode->pointer1)->point[u].y,price};
     // cout<<(treeNode->pointer1)->point[u].x<<" "<<(treeNode->pointer1)->point[u].y<<" "<<price<<"\n";
      finalOutput.push(ans);
    }

    }
    else
    {
     // cout<<"In First";
      calculateMinMaxPrice(xval,yval,treeNode->pointer1);
      q.push(treeNode->pointer1);
     
    }
  }

  if(treeNode->pointer2!=NULL)
  {
    if((treeNode->pointer2)->pointer1==NULL && (treeNode->pointer2)->pointer2==NULL && (treeNode->pointer2)->pointer3==NULL && (treeNode->pointer2)->pointer4==NULL)
    {
      //cout<<"Virti";
     for(int u=0;u<(treeNode->pointer2)->point.size();u++)
      {
      double euclid = sqrt(pow((treeNode->pointer2)->point[u].x-xval,2)+pow((treeNode->pointer2)->point[u].y-yval,2));
      double price = (treeNode->pointer2)->point[u].w - euclid;
      Test ans = {(treeNode->pointer2)->point[u].x,(treeNode->pointer2)->point[u].y,price};
      //cout<<(treeNode->pointer2)->point[u].x<<" "<<(treeNode->pointer2)->point[u].y<<" "<<price<<"\n";
      finalOutput.push(ans);
    }
    }
      else
    {
      calculateMinMaxPrice(xval,yval,treeNode->pointer2);
      q.push(treeNode->pointer2);
      
    }
  }

  if(treeNode->pointer3!=NULL)
  {
    if((treeNode->pointer3)->pointer1==NULL && (treeNode->pointer3)->pointer2==NULL && (treeNode->pointer3)->pointer3==NULL && (treeNode->pointer3)->pointer4==NULL)
    {
      //cout<<"Navil";
      for(int u=0;u<(treeNode->pointer3)->point.size();u++)
      {
      double euclid = sqrt(pow((treeNode->pointer3)->point[u].x-xval,2)+pow((treeNode->pointer3)->point[u].y-yval,2));
      double price = (treeNode->pointer3)->point[u].w - euclid;
      Test ans = {(treeNode->pointer3)->point[u].x,(treeNode->pointer3)->point[u].y,price};
     //cout<<(treeNode->pointer3)->point[u].x<<" "<<(treeNode->pointer3)->point[u].y<<" "<<price<<"\n";
      finalOutput.push(ans);
    }
    }
        else
    {
      calculateMinMaxPrice(xval,yval,treeNode->pointer3);
      q.push(treeNode->pointer3);
      
    
    }
  }
  if(treeNode->pointer4)
  {
    if((treeNode->pointer4)->pointer1==NULL && (treeNode->pointer4)->pointer2==NULL && (treeNode->pointer4)->pointer3==NULL && (treeNode->pointer4)->pointer4==NULL)
    {
       //cout<<"Parita";
      for(int u=0;u<(treeNode->pointer4)->point.size();u++)
      {
      double euclid = sqrt(pow((treeNode->pointer4)->point[u].x-xval,2)+pow((treeNode->pointer4)->point[u].y-yval,2));
      double price = (treeNode->pointer4)->point[u].w - euclid;
      Test ans = {(treeNode->pointer4)->point[u].x,(treeNode->pointer4)->point[u].y,price};
      //cout<<(treeNode->pointer4)->point[u].x<<" "<<(treeNode->pointer4)->point[u].y<<" "<<price<<"\n";
      finalOutput.push(ans);
    }
    }
    else
    {
      calculateMinMaxPrice(xval,yval,treeNode->pointer4);
      q.push(treeNode->pointer4);
    
    }
   }
    int countVal=0;
    vector<double> temp;
    
    while(!q.empty() && countVal<k)
    {
      struct Node* a = q.top();
     // cout<<"In start   ";
      //cout<<"Val = "<<a->xmin<<" "<<a->ymin<<" "<<" "<<a->xmax<<" "<<a->ymax<<" "<<a->count<<"\n";
      temp.push_back(a->minprice);
      q.pop();
      finalQueue.push(a);
      countVal = countVal+a->count;
    }
   //cout<<"Helooooooo\n";
    while(!q.empty())
    {
     // cout<<"HII\n";
       struct Node* b = q.top();
       for(int cin=0;cin<temp.size();cin++)
       {
          if(b->maxprice>=temp[cin])
          {
            //cout<<"Val = "<<b->minprice<<" "<<b->maxprice<<" "<<b->count<<"\n";
            finalQueue.push(b);
            break;
          }
       }

      q.pop();
    }
    //cout<<(q.top())->xmin<<" "<<(q.top())->ymin<<" "<<(q.top())->xmax<<" "<<(q.top())->ymax<<" "<<(q.top())->count<<" "<<(q.top())->minprice<<"\n";
    
    double maxmin = (q.top())->minprice;
     double z;
    while(!finalQueue.empty())
    {
     // cout<<"Gooddyy";
       //cout<<"Final Size = "<<finalQueue.size()<<"\n";
      z = getOutput(finalQueue,maxmin,xval,yval); 
      maxmin = z;
    }
    //cout<<"\n";
   // cout<<"final output size = "<<finalOutput.size()<<"\n";
    while(k!=0)
    {
      Test o = finalOutput.top();
       //cout<<o.x<<" "<<o.y<<"\n";
      myfile2<<o.x<<" "<<o.y<<"\n";
      finalOutput.pop();
      k--;
    }
  auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start1); 
      myfile2<< "\nTime taken by function: "<< duration.count() << " microseconds" ;
      myfile2<<"\n";
      timefile2<<duration.count()<<"\n";
    //  q.clear();
      //finalOutput.clear();
      temp.clear();
      while(!finalOutput.empty())
      {
        finalOutput.pop();
      }
      while(!q.empty())
      {
        q.pop();
      }

}
    




    std::ifstream iqs("queries5.dt.cin");
  std::ofstream myfile5;
  std::ofstream timefile5;
  myfile5.open("FinalOutput5.txt");
  timefile5.open("TimeFile5.txt");
 // assert( ils );
 
  while(iqs>>p1)
  {
    int k=5;
    auto start1 = high_resolution_clock::now(); 
    ios_base::sync_with_stdio(false);
    xval=p1.x();
    yval=p1.y(); 
    myfile5<<"-----------------------------------------------------------------------\n";
    myfile5<<"*****The top-"<<k<<" ads for points"<<xval<<" "<<yval<<" are ****"<<"\n";
    priority_queue<Node*,vector<Node*>,compare> q;

    struct Node* treeNode = head;
    //struct Node* node1 = treeNode->pointer1;
    if(treeNode->pointer1!=NULL)
    {
    if((treeNode->pointer1)->pointer1==NULL && (treeNode->pointer1)->pointer2==NULL && (treeNode->pointer1)->pointer3==NULL && (treeNode->pointer1)->pointer4==NULL)
    {
         for(int u=0;u<(treeNode->pointer1)->point.size();u++)
      {
       double euclid = sqrt(pow((treeNode->pointer1)->point[u].x-xval,2)+pow((treeNode->pointer1)->point[u].y-yval,2));
       double price = (treeNode->pointer1)->point[u].w - euclid;
      Test ans = {(treeNode->pointer1)->point[u].x,(treeNode->pointer1)->point[u].y,price};
      //cout<<(treeNode->pointer1)->point[u].x<<" "<<(treeNode->pointer1)->point[u].y<<" "<<price<<"\n";
      finalOutput.push(ans);
    }

    }
    else
    {
      //cout<<"In First";
      calculateMinMaxPrice(xval,yval,treeNode->pointer1);
      q.push(treeNode->pointer1);
     
    }
  }

  if(treeNode->pointer2!=NULL)
  {
    if((treeNode->pointer2)->pointer1==NULL && (treeNode->pointer2)->pointer2==NULL && (treeNode->pointer2)->pointer3==NULL && (treeNode->pointer2)->pointer4==NULL)
    {
      //cout<<"Virti";
     for(int u=0;u<(treeNode->pointer2)->point.size();u++)
      {
       double euclid = sqrt(pow((treeNode->pointer2)->point[u].x-xval,2)+pow((treeNode->pointer2)->point[u].y-yval,2));
       double price = (treeNode->pointer2)->point[u].w - euclid;
      Test ans = {(treeNode->pointer2)->point[u].x,(treeNode->pointer2)->point[u].y,price};
     // cout<<(treeNode->pointer2)->point[u].x<<" "<<(treeNode->pointer2)->point[u].y<<" "<<price<<"\n";
      finalOutput.push(ans);
    }
    }
      else
    {
      calculateMinMaxPrice(xval,yval,treeNode->pointer2);
      q.push(treeNode->pointer2);
      
    }
  }
  if(treeNode->pointer3!=NULL)
  {
    if((treeNode->pointer3)->pointer1==NULL && (treeNode->pointer3)->pointer2==NULL && (treeNode->pointer3)->pointer3==NULL && (treeNode->pointer3)->pointer4==NULL)
    {
      //cout<<"Navil";
      for(int u=0;u<(treeNode->pointer3)->point.size();u++)
      {
       double euclid = sqrt(pow((treeNode->pointer3)->point[u].x-xval,2)+pow((treeNode->pointer3)->point[u].y-yval,2));
       double price = (treeNode->pointer3)->point[u].w - euclid;
      Test ans = {(treeNode->pointer3)->point[u].x,(treeNode->pointer3)->point[u].y,price};
     // cout<<(treeNode->pointer3)->point[u].x<<" "<<(treeNode->pointer3)->point[u].y<<" "<<price<<"\n";
      finalOutput.push(ans);
    }
    }
        else
    {
      calculateMinMaxPrice(xval,yval,treeNode->pointer3);
      q.push(treeNode->pointer3);
      
    
    }
  }
  if(treeNode->pointer4!=NULL)
  {
    if((treeNode->pointer4)->pointer1==NULL && (treeNode->pointer4)->pointer2==NULL && (treeNode->pointer4)->pointer3==NULL && (treeNode->pointer4)->pointer4==NULL)
    {
      //  cout<<"Parita";
      for(int u=0;u<(treeNode->pointer4)->point.size();u++)
      {
       double euclid = sqrt(pow((treeNode->pointer4)->point[u].x-xval,2)+pow((treeNode->pointer4)->point[u].y-yval,2));
       double price = (treeNode->pointer4)->point[u].w - euclid;
      Test ans = {(treeNode->pointer4)->point[u].x,(treeNode->pointer4)->point[u].y,price};
      //cout<<(treeNode->pointer4)->point[u].x<<" "<<(treeNode->pointer4)->point[u].y<<" "<<price<<"\n";
      finalOutput.push(ans);
    }
    }
    else
    {
      calculateMinMaxPrice(xval,yval,treeNode->pointer4);
      q.push(treeNode->pointer4);
    
    }
   }
    
    int countVal=0;
    vector< double> temp;
    
    while(!q.empty() && countVal<k)
    {
      struct Node* a = q.top();
     // cout<<"In start   ";
      //cout<<"Val = "<<a->xmin<<" "<<a->ymin<<" "<<" "<<a->xmax<<" "<<a->ymax<<" "<<a->count<<"\n";
      temp.push_back(a->minprice);
      q.pop();
      finalQueue.push(a);
      countVal = countVal+a->count;
    }
   //cout<<"Helooooooo\n";
    while(!q.empty())
    {
     // cout<<"HII\n";
       struct Node* b = q.top();
       for(int cin=0;cin<temp.size();cin++)
       {
          if(b->maxprice>=temp[cin])
          {
            //cout<<"Val = "<<b->minprice<<" "<<b->maxprice<<" "<<b->count<<"\n";
            finalQueue.push(b);
            break;
          }
       }

      q.pop();
    }
    //cout<<(q.top())->xmin<<" "<<(q.top())->ymin<<" "<<(q.top())->xmax<<" "<<(q.top())->ymax<<" "<<(q.top())->count<<" "<<(q.top())->minprice<<"\n";
    
     double maxmin = (q.top())->minprice;
     double z;
    while(!finalQueue.empty())
    {
     // cout<<"Gooddyy";
       //cout<<"Final Size = "<<finalQueue.size()<<"\n";
      z = getOutput(finalQueue,maxmin,xval,yval); 
      maxmin = z;
    }
    //cout<<"\n";
    //cout<<"final output size = "<<finalOutput.size()<<"\n";
    while(k!=0)
    {
      Test o = finalOutput.top();
      myfile5<<o.x<<" "<<o.y<<"\n";
      finalOutput.pop();
      k--;
    }
  auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start1); 
      myfile5<< "\nTime taken by function: "<< duration.count() << " microseconds" ;
      myfile5<<"\n";
      timefile5<<duration.count()<<"\n";
    //  q.clear();
      //finalOutput.clear();
      temp.clear();
      while(!finalOutput.empty())
      {
        finalOutput.pop();
      }
      while(!q.empty())
      {
        q.pop();
      }


    }


 std::ifstream ips("queries10.dt.cin");
  std::ofstream myfile10;
  std::ofstream timefile10;
  myfile10.open("FinalOutput10.txt");
  timefile10.open("TimeFile10.txt");
 // assert( ils );
 
  while(ips>>p1)
  {
    int k=10;
    auto start1 = high_resolution_clock::now(); 
    ios_base::sync_with_stdio(false);
    xval=p1.x();
    yval=p1.y(); 
    myfile10<<"-----------------------------------------------------------------------\n";
    myfile10<<"*****The top-"<<k<<" ads for points"<<xval<<" "<<yval<<" are ****"<<"\n";
    priority_queue<Node*,vector<Node*>,compare> q;

    struct Node* treeNode = head;
    //struct Node* node1 = treeNode->pointer1;
    if(treeNode->pointer1!=NULL)
    {
    if((treeNode->pointer1)->pointer1==NULL && (treeNode->pointer1)->pointer2==NULL && (treeNode->pointer1)->pointer3==NULL && (treeNode->pointer1)->pointer4==NULL)
    {
         for(int u=0;u<(treeNode->pointer1)->point.size();u++)
      {
       double euclid = sqrt(pow((treeNode->pointer1)->point[u].x-xval,2)+pow((treeNode->pointer1)->point[u].y-yval,2));
       double price = (treeNode->pointer1)->point[u].w - euclid;
      Test ans = {(treeNode->pointer1)->point[u].x,(treeNode->pointer1)->point[u].y,price};
      //cout<<(treeNode->pointer1)->point[u].x<<" "<<(treeNode->pointer1)->point[u].y<<" "<<price<<"\n";
      finalOutput.push(ans);
    }

    }
    else
    {
      //cout<<"In First";
      calculateMinMaxPrice(xval,yval,treeNode->pointer1);
      q.push(treeNode->pointer1);
     
    }
    }

    if(treeNode->pointer2!=NULL)
    {
    if((treeNode->pointer2)->pointer1==NULL && (treeNode->pointer2)->pointer2==NULL && (treeNode->pointer2)->pointer3==NULL && (treeNode->pointer2)->pointer4==NULL)
    {
      //cout<<"Virti";
     for(int u=0;u<(treeNode->pointer2)->point.size();u++)
      {
       double euclid = sqrt(pow((treeNode->pointer2)->point[u].x-xval,2)+pow((treeNode->pointer2)->point[u].y-yval,2));
       double price = (treeNode->pointer2)->point[u].w - euclid;
      Test ans = {(treeNode->pointer2)->point[u].x,(treeNode->pointer2)->point[u].y,price};
     // cout<<(treeNode->pointer2)->point[u].x<<" "<<(treeNode->pointer2)->point[u].y<<" "<<price<<"\n";
      finalOutput.push(ans);
    }
    }
      else
    {
      calculateMinMaxPrice(xval,yval,treeNode->pointer2);
      q.push(treeNode->pointer2);
      
    }
    }

    if(treeNode->pointer3!=NULL)
    {
    if((treeNode->pointer3)->pointer1==NULL && (treeNode->pointer3)->pointer2==NULL && (treeNode->pointer3)->pointer3==NULL && (treeNode->pointer3)->pointer4==NULL)
    {
      //cout<<"Navil";
      for(int u=0;u<(treeNode->pointer3)->point.size();u++)
      {
       double euclid = sqrt(pow((treeNode->pointer3)->point[u].x-xval,2)+pow((treeNode->pointer3)->point[u].y-yval,2));
       double price = (treeNode->pointer3)->point[u].w - euclid;
      Test ans = {(treeNode->pointer3)->point[u].x,(treeNode->pointer3)->point[u].y,price};
     // cout<<(treeNode->pointer3)->point[u].x<<" "<<(treeNode->pointer3)->point[u].y<<" "<<price<<"\n";
      finalOutput.push(ans);
    }
    }
        else
    {
      calculateMinMaxPrice(xval,yval,treeNode->pointer3);
      q.push(treeNode->pointer3);
      
    
    }
  }

  if(treeNode->pointer4!=NULL)
  {
    if((treeNode->pointer4)->pointer1==NULL && (treeNode->pointer4)->pointer2==NULL && (treeNode->pointer4)->pointer3==NULL && (treeNode->pointer4)->pointer4==NULL)
    {
      //  cout<<"Parita";
      for(int u=0;u<(treeNode->pointer4)->point.size();u++)
      {
       double euclid = sqrt(pow((treeNode->pointer4)->point[u].x-xval,2)+pow((treeNode->pointer4)->point[u].y-yval,2));
       double price = (treeNode->pointer4)->point[u].w - euclid;
      Test ans = {(treeNode->pointer4)->point[u].x,(treeNode->pointer4)->point[u].y,price};
      //cout<<(treeNode->pointer4)->point[u].x<<" "<<(treeNode->pointer4)->point[u].y<<" "<<price<<"\n";
      finalOutput.push(ans);
    }
    }
    else
    {
      calculateMinMaxPrice(xval,yval,treeNode->pointer4);
      q.push(treeNode->pointer4);
    
    }
   
    }
    int countVal=0;
    vector< double> temp;
    
    while(!q.empty() && countVal<k)
    {
      struct Node* a = q.top();
     // cout<<"In start   ";
      //cout<<"Val = "<<a->xmin<<" "<<a->ymin<<" "<<" "<<a->xmax<<" "<<a->ymax<<" "<<a->count<<"\n";
      temp.push_back(a->minprice);
      q.pop();
      finalQueue.push(a);
      countVal = countVal+a->count;
    }
   //cout<<"Helooooooo\n";
    while(!q.empty())
    {
     // cout<<"HII\n";
       struct Node* b = q.top();
       for(int cin=0;cin<temp.size();cin++)
       {
          if(b->maxprice>=temp[cin])
          {
            //cout<<"Val = "<<b->minprice<<" "<<b->maxprice<<" "<<b->count<<"\n";
            finalQueue.push(b);
            break;
          }
       }

      q.pop();
    }
    //cout<<(q.top())->xmin<<" "<<(q.top())->ymin<<" "<<(q.top())->xmax<<" "<<(q.top())->ymax<<" "<<(q.top())->count<<" "<<(q.top())->minprice<<"\n";
    
     double maxmin = (q.top())->minprice;
     double z;
    while(!finalQueue.empty())
    {
     // cout<<"Gooddyy";
       //cout<<"Final Size = "<<finalQueue.size()<<"\n";
      z = getOutput(finalQueue,maxmin,xval,yval); 
      maxmin = z;
    }
    //cout<<"\n";
    //cout<<"final output size = "<<finalOutput.size()<<"\n";
    while(k!=0)
    {
      Test o = finalOutput.top();
      myfile10<<o.x<<" "<<o.y<<"\n";
      finalOutput.pop();
      k--;
    }
  auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start1); 
      myfile10<< "\nTime taken by function: "<< duration.count() << " microseconds" ;
      myfile10<<"\n";
      timefile10<<duration.count()<<"\n";
    //  q.clear();
      //finalOutput.clear();
      temp.clear();
      while(!finalOutput.empty())
      {
        finalOutput.pop();
      }
      while(!q.empty())
      {
        q.pop();
      }


    }




	return 0;
}