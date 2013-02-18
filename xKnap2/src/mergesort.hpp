// QUELLE: http://ankitstar.blogspot.de/2011/05/merge-sort.html

void merge(vector<double>*,vector<int>*,int,int,int);

void mergesort(vector<double>* a,vector<int>* merke,int p,int r) // a=Vector of length n ,p=0,r=n-1
{
     if( p < r)
    {
         int q=(p+r)/2;
         mergesort(a,merke,p,q);
         mergesort(a,merke,q+1,r) ;
         merge(a,merke,p,q,r);
     }
}
void merge(vector<double>* a,vector<int>* merke,int p,int q,int r)
{
double c[10];
double d[10];
int i=p;
int j=q+1;
int k=p;
while((i<=q)&&(j<=r))
{
  if((*a).at(i)<(*a).at(j))
  {
    c[k]=(*a).at(i);
    d[k]=(*merke).at(i);
    i=i+1;
    k=k+1;
  }
  else
  {
    c[k]=(*a).at(j);
    d[k]=(*merke).at(j);
    j=j+1;
    k=k+1;
  }
}
while(i<=q)
{
  c[k] =(*a).at(i);
  d[k]=(*merke).at(i);
  i=i+1;
  k=k+1;
}
while(j<=r)
{
  c[k]=(*a).at(j);
  d[k]=(*merke).at(j);
  j=j+1;
  k=k+1;
}
int l=p;
while(l<=r)
{
  (*a).at(l)=c[l];
  (*merke).at(l)=d[l];
  l=l+1;
}
}
