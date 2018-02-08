//Khalid & Vidya

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>

/*******************************************************************************************/
//double* sort_cluster;
double* get_random_datapoint(double* data, int ndata, int dim, int cluster_id, int* cluster_size,int* cluster_start,int *cluster_assign);
double* get_next_centroid(double* data,int dim, double* point, int cluster_id, int* cluster_assign, int ndata,int *cluster_size);
double get_mean_of_dim(double* data,int *cluster_assign, int* cluster_size, int current_cluster_index,int dim_id, int ndata);
double* find_centroid_for_empty_cluster(double* data, int* cluster_assign,int cluster_id, double* point, int ndata, int dim );
double calculate_eulidean_distance(double* point1, double* point2, int dim);
/*********************************************************************************************/
void bisecting_kmeans(int dim, int ndata,
                      double *data,
                      int k,
                      int *cluster_size,
                      int *cluster_start,
                      double *cluster_radius,
                      double *cluster_centroid,
                      int *cluster_assign) {
    int i;
    printf("-----------------------Bisecting initial cluster--------------------------------------------  \n");
    int current_number_of_clusters = 1, cluster_id_of_max_SSE = 0, size = 0;
    printf("-------------------------Printing data and cluster assign-----------------------------------\n");
    print_data(data, cluster_assign, ndata, dim);


    while (current_number_of_clusters < k) //iterative process of bisecting kmeans until we get k no of clusters.
    {
        if (cluster_size[cluster_id_of_max_SSE] != 0) //call of bisecting kmeans functions to the every first cluster.
        {
            double *first_centroid = get_random_datapoint(data, ndata, dim, cluster_id_of_max_SSE, cluster_size,
                                                          cluster_start,
                                                          cluster_assign);//get random datapoint for the cluster as a initial centroid
            double *second_centroid = get_next_centroid(data, dim, first_centroid, cluster_id_of_max_SSE,
                                                        cluster_assign, ndata,
                                                        cluster_size); //get next centroid of the cluster by calculating the farthest distance from initial centroid
            //assign each datapoint to nearest cluster
            nearest_cluster_assign(first_centroid, second_centroid, data, ndata, dim, cluster_assign, &cluster_assign,
                                   cluster_size, cluster_centroid, cluster_id_of_max_SSE, current_number_of_clusters, k,
                                   &cluster_centroid, &cluster_size);

            //calculate the maximum SSE and retrieve the cluster id of max SSE.
            cluster_id_of_max_SSE = get_cluster_id_of_max_SSE(cluster_assign, data, cluster_id_of_max_SSE,
                                                              current_number_of_clusters, cluster_centroid, dim, ndata);

        }
        //increment the no of clusters
        current_number_of_clusters++;

    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    printf("==================================================After bisecting is done======================================\n");
    //sort the data according to cluster assign finally when we get k no of clusters as requested.
    printf("----------------------------------------------sorting of data points------------------------------------\n");
    sort_cluster(cluster_assign, &cluster_assign, data, &data, &cluster_start, ndata, dim, k,
                 cluster_size, &cluster_size);
    printf("Printing data and cluster assign after sorting the data only:\n");
    print_data(data, cluster_assign, ndata, dim); //print the sorted data points.
    for(i=0; i<ndata; i++)
        printf("cluster_assign[%i]= %i\n",i,cluster_assign[i]);
    printf("-------------------------------cluster_start marks-------------------------------------------------------\n");
    print_normarl_array(cluster_start, k, 0); //print the start of every cluster formed.
    printf("---------------------------------cluster_size recorded---------------------------------------------------\n");
    print_normarl_array(cluster_size, k, 1); //prints the cluster_size of every cluster formed.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (i = 0; i < k * dim; i++)
        printf("cluster_centroid[%i]= %lf\n", i, cluster_centroid[i]); //print the centroid of all k clusters.

///////////////////////////////////////////////////////////////////////////////////////////////////
    //function call to calculate the radius of each clusters.
    find_cluster_radius(data, cluster_size, cluster_start, ndata, dim, k, cluster_centroid, &cluster_radius);

    for(i=0;i<k;i++)
        printf("cluster_radius[%i]: %lf\n",i,cluster_radius[i]);

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    //function to search the closest point to a randomly generated query point.
    //  Search(dim, ndata, data, k, cluster_radius, cluster_centroid, cluster_assign, cluster_size, cluster_start);


}
/********************************************************************************************************/
/*calculate the mean of every dimension of a data point in order calculate the SSE and to check the static cluster assign*/
double get_mean_of_dim(double* data,int *cluster_assign, int* cluster_size, int current_cluster_index,int dim_id, int ndata)
{
    int i;
    double sum = 0;
    ////////////////////////////////////////////////////////////////////////////////////
    if(cluster_size[current_cluster_index]!=0) {
        for (i = 0; i < ndata; i++) {
            if (cluster_assign[i] == current_cluster_index) {
                //sum the data points according to the dimension value that we want mean for
                sum = sum + data[i + dim_id];
            }
        }
        //divide the sum by no of clusters
        return (sum / cluster_size[current_cluster_index]);
    }
    else{
        //if cluster size is empty then mean is randomly generated.
        return (rand()%ndata-1);
    }
}
/********************************************************************************************************/
/*function to find the nearest cluster*/
int nearest_cluster_assign(double* first_centroid,double* next_centroid,double* data,int ndata,int dim,int* cluster_assign,int** out_cluster_assign,int* cluster_size,double* cluster_centroid, int cluster_id, int no_of_clusters_currently,int k,double **out_cluster_centroid, int **out_cluster_size){


    int count=1,count1=1, flag=0, counter=0, index=0;
    double dist1,dist2;
    int *temp_cluster_assign = cluster_assign;
    int *temp_cluster_assign2 = cluster_assign;
    int *temp_cluster_size = cluster_size;
    double * temp_cluster_centroid= cluster_centroid;
    int i,j;
    int enterance_flag_first_cluster=0, enterance_flag_second_cluster=0;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    while(flag==0)
    {
       // if(cluster_size[cluster_id]!=0)
        {
            for (i = 0; i < ndata; i++) {
                //for the particular cluster
                if (cluster_assign[i] == cluster_id) {
                    //calculate the distance from both the first and the next centroid to every datapoint
                    dist1 = calculate_eulidean_distance(first_centroid, &data[get_location_of_datapoint(i, dim)], dim);

                    dist2 = calculate_eulidean_distance(next_centroid, &data[get_location_of_datapoint(i, dim)], dim);
                    //compare the distance values
                    if (dist1 < dist2) {
                        enterance_flag_first_cluster=1;
                        //if dist1 is less than dist2 then assign the datapoint to the first centroids cluster
                        temp_cluster_assign[i] = cluster_id;
                        temp_cluster_size[cluster_id] = count;  //increase the cluster size by 1
                        count++;
                    } else {

                        enterance_flag_second_cluster=1;
                        //if dist2 is else then dist1 then assign the datapoint to the next centroids cluster
                        temp_cluster_assign[i] = no_of_clusters_currently;
                        temp_cluster_size[no_of_clusters_currently] = count1; //increment the cluster size by 1
                        count1++;
                    }
                }
            }
            //////////////////////////updating centroids
            counter = 0;
            /*the entrance flags help us detecting the empty cluster instead of looping through the cluste size array everytime
         * initially those flags are zeros, but if we started filling the cluster size array for that cluste on the above code, the flags will be one*/
            //if both flags = 1, that means both clusters have data and are not empty
            if(enterance_flag_first_cluster==1 && enterance_flag_second_cluster==1) {
                for (j = 0; j < dim; j++) {

                    index = (cluster_id * dim) + j;
                    temp_cluster_centroid[index] = get_mean_of_dim(data, temp_cluster_assign, temp_cluster_size,
                                                                   cluster_id, j, ndata);
                    first_centroid[j] = temp_cluster_centroid[index];
                }


                // get the centroid of second cluster
                for (j = 0; j < dim; j++) {
                    index = (no_of_clusters_currently * dim) + j;
                    temp_cluster_centroid[index] = get_mean_of_dim(data, temp_cluster_assign, temp_cluster_size,
                                                                   no_of_clusters_currently, j, ndata);
                    next_centroid[j] = temp_cluster_centroid[index];
                }


            }
                // if the flag of the first cluster is one, it means the other cluster is empty
            else if(enterance_flag_first_cluster==1 && enterance_flag_second_cluster==0)
            {
                //in this case we calculate the mean of the first cluster only and get the centroid
                for (j = 0; j < dim; j++) {

                    index = (cluster_id * dim) + j;
                    temp_cluster_centroid[index] = get_mean_of_dim(data, temp_cluster_assign, temp_cluster_size,
                                                                   cluster_id, j, ndata);
                    first_centroid[j] = temp_cluster_centroid[index];
                }

                //then we find a mean for the other cluster using the previous centroid
                next_centroid = find_centroid_for_empty_cluster( data,  cluster_assign,cluster_id,  first_centroid,  ndata,dim );
                for (j = 0; j < dim; j++) {
                    index = (no_of_clusters_currently * dim) + j;
                    temp_cluster_centroid[index] = next_centroid[j];

                }
            }
                // if the flag of the second cluster is one, it means the first cluster is empty
            else if(enterance_flag_first_cluster==0 && enterance_flag_second_cluster==1)
            {
                //in this case we calculate the mean of the second centroid first because the original cluster became empty after bisecting
                for (j = 0; j < dim; j++) {
                    index = (no_of_clusters_currently * dim) + j;
                    temp_cluster_centroid[index] = get_mean_of_dim(data, temp_cluster_assign, temp_cluster_size,
                                                                   no_of_clusters_currently, j, ndata);
                    next_centroid[j] = temp_cluster_centroid[index];
                }
                //then we find a mean for the other cluster using the previous centroid
                first_centroid = find_centroid_for_empty_cluster( data,cluster_assign,no_of_clusters_currently,next_centroid, ndata,dim );
                for (j = 0; j < dim; j++) {
                    index = (cluster_id * dim) + j;
                    temp_cluster_centroid[index] = first_centroid[j];
                }
            }
            flag = check_validity_of_cluster_assign(temp_cluster_assign2, temp_cluster_assign, ndata);
            temp_cluster_assign2 = temp_cluster_assign;
        }
        //else
        {
          //  flag=1;
        }
    }

    *out_cluster_assign=temp_cluster_assign2;
    *out_cluster_centroid=temp_cluster_centroid;
    *out_cluster_size=temp_cluster_size;
}
/********************************************************************************************************/
/* function to find the centroid of the empty cluster*/
double* find_centroid_for_empty_cluster(double* data, int* cluster_assign,int cluster_id, double* point, int ndata, int dim )
{
    int nextcentroid_index,i;
    double next_max_dist=0, max_dist=0;
    double *result_point = NULL;
    /////////////////////////////////////////////////////////////////////
    for(i=0;i<ndata;i++)
    {
        if (cluster_assign[i] == cluster_id) {
            //calculate the distance from current centroid to all the other datapoints and the one datapoint with maximum
            //distance will be the next centroid
            next_max_dist = calculate_eulidean_distance(point, &data[get_location_of_datapoint(i, dim)], dim);
            if (next_max_dist >= max_dist) {
                max_dist = next_max_dist;//record the maximum distance found
                nextcentroid_index = i; //record the datapoint of maximum distance
                //now this datapoint will be the centroid of the empty cluster
            }
        }
    }
    //store the datapoint into pointer
    result_point = &data[get_location_of_datapoint(nextcentroid_index,dim)];
    //return the new centroid for empty cluster
    return result_point;
}
/********************************************************************************************************/
/*function to calculate the max SSE of a cluster and also to return the cluster_id with max SSE*/
int get_cluster_id_of_max_SSE(int *cluster_assign, double* data, int cluster_1, int cluster_2, double* cluster_centroid, int dim, int ndata)
{
    double SSE_cluster1 = 0, SSE_cluster2=0;
    int i,j;
    //////////////////////////////////////////////////////////////////
    for( i=0; i<ndata;i++)
    {
        for( j=0; j<dim;j++)
        {
            //each time the SSE is calculated between two clusters that are needed next to bisect
            if(cluster_assign[i]==cluster_1)
                //calculate the SSE based on the formula
                SSE_cluster1=SSE_cluster1+sqrt(data[i+j]- cluster_centroid[cluster_1+j]);
            if(cluster_assign[i]==cluster_2)
                SSE_cluster2=SSE_cluster2+sqrt(data[i+j]- cluster_centroid[cluster_2+j]);
        }
    }
    /////////////////////////////////////////////////////////////////////////
    //compare the SSE and return the maximum SSE cluster no
    if(SSE_cluster1>SSE_cluster2)
        return cluster_1;
    else
        return cluster_2;
}
/********************************************************************************************************/
/* function to check the validity of cluster*/
//to check whether the data point belongs to the same cluster or not
int check_validity_of_cluster_assign(int * cluster_assign, int * temp_cluster_assign, int ndata)
{
    int i;
    for( i=0; i<ndata; i++)
    {
        if(cluster_assign[i]!=temp_cluster_assign[i])
            return 0;
    }
    return 1;
}
/********************************************************************************************************/
/*function to calculate the cluster radius*/
int find_cluster_radius(double* data,int* cluster_size,int* cluster_start,int ndata,int dim,int k,double* cluster_centroid,double** cluster_radius){
    int i,j,start;
    double * max_radius = (double *)malloc(sizeof(double) * k* dim);
    /////////////////////////////////////////////////////////////////////
    for ( i = 0; i < k; i++) {  //array to store the cluster radius
        max_radius[i] = 0;
    };
    ////////////////////////////////////////////////////////////////////
    printf("===================================cluster radius=======================================\n");
    double distance = 0;
    //cluster radius is the distance between the cluster centroid and the furthest datapoint
    for(i=0;i<k;i++){
        for(j=0;j<(cluster_size[i]);j++){ //loop till the cluster size
            start=cluster_start[i]; //begin from the cluster start of the particular cluster that you want to find radius
            //printf("start+j:%d\n",start+j);
            //cluster centroid is given to pointer point
            double* point = &cluster_centroid[get_location_of_datapoint(i,dim)];
            //data is assigned to pointer point1
            double* point1 = &data[start+j*dim];
            //calculate the distance from the centroid to each data of a cluster
            distance = calculate_eulidean_distance(point,point1,dim);
            // printf("radius distance:%lf\n",distance);
            if (distance > max_radius[i])
            {
                max_radius[i] = distance; //record the maximum distance found
            }
        }
    }
    /////////////////////////////////////////////////////////////////////////
    // for (i = 0; i < k; i++)
    //{
    *cluster_radius= max_radius; //store the maximum distance in cluster radius
    //}
    /////////////////////////////////////////////////////////////////////

}

/**********************************************************************************************************/
/*sort the data based on the cluster_assign and also retrieve the cluster start and size accordingly*/
int sort_cluster(int* cluster_assign, int** out_cluster_assign, double *data,double **dataOut, int **cluster_start, int ndata, int dim, int k,int *cluster_size, int **out_cluster_size){
    double *temp_Data= (double *)malloc(sizeof(double) * ndata * dim); //temp data set array
    int *temp_cluster_assign= (int *)malloc(sizeof(int) * ndata);//temp cluster assign
    int *temp_cluster_start= (int *)malloc(sizeof(int) * k); //temp cluster start
    int *temp_cluster_size= (int *)malloc(sizeof(int) * k);
    int data_index=0, assign_index=0, start_flag=0, size_counter=0;
    int cluster,i,j;
    ///////////////////////////////////////////////////////////////////////////////////
    for( cluster=0; cluster<k;cluster++) { //k no of clusters created
        start_flag = 0;

        for (i = 0; i < ndata; i++) { //for ndata of datapoints
            if (cluster_assign[i] == cluster) { //cluster assign will be the cluster no
                size_counter++;  //for every datapoint of a cluster the size counter increments
                if (start_flag == 0) {
                    if (cluster_size[cluster] != 0)
                        temp_cluster_start[cluster] = data_index; //if cluster size is !=0 then the cluster start will be the dataindex of the
                    //first datapoint of that cluster
                    start_flag = 1;
                }

                for (j = 0; j < dim; j++) {
                    //store the data into the temp data array
                    temp_Data[data_index + j] = data[get_location_of_datapoint(i, dim) + j];
                }
                data_index = data_index + dim;
                //store the cluster assign or cluster no to the temp cluster assign
                temp_cluster_assign[assign_index++] = cluster;
            }
        }
        //store the cluster size into temp cluster size
        temp_cluster_size[cluster] = size_counter;
        size_counter = 0;
    }
    //////////////////////////////////////////////////////////////
    //return all the temp array of data,cluster size,cluster assign and cluster start/
    *dataOut=temp_Data;
    *cluster_start=temp_cluster_start;
    *out_cluster_assign=temp_cluster_assign;
    *out_cluster_size=temp_cluster_size;
}

/****************************************************************************************************/
//print the cluster size,cluster start and cluster centroids
int  print_normarl_array(int * data,int k, int flag)
{
    int j;
    for(  j=0; j<k; j++)
    {
        if(flag==0)
            printf("cluster_start[%i]: %i\n",j,data[j]);
        if(flag==1)
            printf("cluster_size[%i]: %i\n",j,data[j]);
        if(flag==2)
            printf("cluster_centroids: %lf\n",(double)data[j]);

    }

}
/****************************************************************************************************/
/*function to print the datapoint and cluster number */
int  print_data(double* data,int* cluster_assign,int ndata,int dim)
{
    int i,j;
    ///////////////////////////////////////////////////
    //print the data point
    for(  i=0; i<ndata; i++)                                          //this for loop print all the datapoints
    {
        printf("(");
        for(  j=0; j<dim; j++)
        {
            printf("%lf",data[get_location_of_datapoint(i,dim)+j]);
            if(j<dim-1)
                printf(", ");
        }
        printf(") \t");
        ///////////////////////////////////////////////////////////////
        //print the cluster no of a particular datapoint
        printf("cluster number (%i)\n",cluster_assign[i]);
    }
}
/*******************************************************************************************/
/*this function return a random datapoint pointer out of our dataset
 * it used in the initial stage when we create the bisecting kmeans*/
double* get_random_datapoint(double* data, int ndata, int dim, int cluster_id, int* cluster_size,int* cluster_start,int *cluster_assign)
{
    int random_index,i;
    double *point = NULL;
    ///////////////////////////////////////////////////////////////////////////
    if(cluster_size[cluster_id]!=0)  //enter if loop if cluster size is !=0
    {
        for( i=0 ; i<ndata; i++)
        {
            random_index = (rand()%ndata-1); //this calculation makes sure not to exceed the the current dataset
            //retrieve the datapoint within a particular cluster with help of random index

            if(cluster_assign[random_index]==cluster_id)
                i=ndata;
        }
        point = &data[get_location_of_datapoint(random_index,dim)];
    }
    else
        //enter else loop if cluster size is 0
    {
        for( i=0; i<dim; i++)
            point[i]= (rand()%ndata-1); //generate the random datapoint as a first centroid
    }
    return point; //return the datapoint as a first centroid
}
/*******************************************************************************************/

/* function to find the next centroid*/
double* get_next_centroid(double* data,int dim, double* point, int cluster_id, int* cluster_assign, int ndata, int* cluster_size)
{
    int nextcentroid_index = 0,i,j;
    double next_max_dist=0, max_dist=0;
    double *result_point = NULL;
    if(cluster_size[cluster_id]!=0) //if cluster size is != 0 that means if cluster size is not empty then enter if loop
    {
        for (i = 0; i < ndata; i++) {
            if (cluster_assign[i] == cluster_id) { //cluster_assign should be cluster_id of maximum SSE that we are working
                //calculate distance from current centroid and all the datapoints within the cluster
                next_max_dist = calculate_eulidean_distance(point, &data[get_location_of_datapoint(i, dim)], dim);
                if (next_max_dist > max_dist) {
                    max_dist = next_max_dist; //record the maximum distance
                    nextcentroid_index = i; //the datapoint with maximum distance is now the next centroid
                }
            }
        }
        result_point = &data[get_location_of_datapoint(nextcentroid_index,dim)]; //return the datapoint as a result point
    }
    else
    {
        for(i=0;i<ndata;i++){
            if (cluster_assign[i] == cluster_id) { //cluster_assign should be cluster_id of maximum SSE that we are working
                //calculate distance from current centroid and all the datapoints within the cluster
                next_max_dist = calculate_eulidean_distance(point, &data[get_location_of_datapoint(i, dim)], dim);
                if (next_max_dist > max_dist) {
                    max_dist = next_max_dist; //record the maximum distance
                    nextcentroid_index = i; //the datapoint with maximum distance is now the next centroid
                }
            }
        }

        result_point = &data[get_location_of_datapoint(nextcentroid_index,dim)]; //return the datapoint as a result point
    }
    return result_point;
}

/****************************************************************************************************/

/* distance between two points */
double calculate_eulidean_distance(double* point1, double* point2, int dim)
{
    double dist = 0;
    int i=0;
    //calculate distance based on euclidean formula of calculating the distance between two points
    for( i = 0; i < dim; i++) {

        dist += (point1[i] - point2[i]) * (point1[i] - point2[i]);
    }
    double distance = sqrt(dist);
    return distance;
}
/*******************************************************************************************/

/* random data generation */
void generate_random_data(double* data, int ndata, int dim) {
    srand(time(NULL));   // should only be called once
    int i;
    for( i = 0; i < (ndata*dim); i++)
    {
        data[i] = (double)((rand() % 100)+1);//random number generator in range of 100
    }
}
/*******************************************************************************************/

/* return the exact location of the element in the array */
int get_location_of_datapoint( int index,int dim) {
    int location=NULL;
    if(index>=0)
    {
        location=dim*index;   // the location holds the index of the datapoint(every dimension)
    }
    return location;
}
/*******************************************************************************************/
/* function to search the nearest datapoint to a randomly generated query point*/
int Search(int dim,int ndata,double *data,int k,double *cluster_radius,double * cluster_centroid, int *cluster_assign,int *cluster_size,int *cluster_start)
{   int i,j;
    double mindist;
    double min_dist=999999999;
    int search_cluster, search_node = 0;
    int start;
    int size;
    double final_dist;
    printf("=========================================search function==================================================\n");
    //////////////////////////////////////////////////////////////////////////
    double* search_query =(double *)malloc(sizeof(double)*dim);

    for(i=0;i<dim;i++){
        search_query[i] =(double)(rand() % 100);
    }
    //print the query point that is randomly generated
    printf("Search node is ");
    printf("(");
    for(i=0;i<dim;i++)
    {
        printf("%lf ",search_query[i]);
        if(j<dim-1)
            printf(",");
    }
    printf(")\n");
    //////////////////////////////////////////////////////////////////////
    //calculate the distance from query point to all the cluster centroids
    for(i=0;i<k;i++){
        if(cluster_size[i]!=0)
        {
            double* point = &cluster_centroid[get_location_of_datapoint(i,dim)];
            mindist = calculate_eulidean_distance(point,&search_query, dim);
            mindist=mindist- cluster_radius[i];
            mindist=abs(mindist);
            if(mindist< min_dist)
            {
                min_dist=mindist;//record the minimum distance occured
                min_dist=abs(min_dist);
                search_cluster=i;//record the cluster no to which we found minimum distance
            }
        }
    }
    printf("The cluster that need to be searched is:%d\n",search_cluster);



    ////////////////////////////////////////////////////////////////////////////
    start=cluster_start[1];
    size=cluster_size[1];
    min_dist=99999999;
    for(i=start ;i< (start+size) ;i+=dim)
    {

        double* point1 = &data[i];
        /* calculate the distance from query point to all the datapoints within the search cluster*/
        final_dist=calculate_eulidean_distance(point1,&search_query,dim);
        // printf("distance is %lf\n",final_dist);
        if(final_dist < min_dist){
            min_dist=final_dist;//record the minimum distance found
            search_node=i;//record the datapoint that is nearer
        }
    }
    //print the datapoint that found nearer to query point
    printf("nearest point found is ");
    printf("(");
    for(i=0;i<dim;i++){
        printf("%lf",data[search_node+i]);
        if(i<dim-1)
            printf(",");
    }
    printf(") \n");
    printf("distance to nearest node is %lf",min_dist);
}
/**********************************************************************************************************/
/*main function*/
int main() {
    int ndata, dim, k, i, j;   // ndata = # of datapoints, dim = # of dimensions k= # of clusters
    /////////////////////////////////////////////////////////////////////////////////////
    printf("Enter the no of Datapoints\n");
    scanf("%d", &ndata);
    printf("Enter the no of Dimension\n");
    scanf("%d", &dim);
    printf("Enter the no of clusters\n");
    scanf("%d", &k);
    /////////////////////////////////////////////////////////////////////////////////////
    /* initializing all the required objects*/
    double *data = (double *) malloc(
            sizeof(double) * ndata * dim);  // define and reserve array for datapoints of double type
    generate_random_data(data, ndata, dim);// generate random datapoints
    //////////////////////////////////////
    int *cluster_size = (int *) malloc(sizeof(int) * k);
    for (i = 0; i < k; i++) {
        if (i == 0)
            cluster_size[i] = ndata;
        else
            cluster_size[i] = 0;
    }
    //////////////////////////////////////
    int *cluster_start = (int *) malloc(sizeof(int) * k);
    for (i = 0; i < k; i++) {
        cluster_start[i] = 0;
    }
    //////////////////////////////////////
    double *cluster_radius = (double *) malloc(sizeof(double) * k - 1 * dim);
    for (i = 0; i < k; i++) {
        cluster_radius[i] = 0;
    }
    //////////////////////////////////////
    double *cluster_centroid = (double *) malloc(sizeof(double) * k * dim);
    for (i = 0; i < k * dim; i++) {
        cluster_centroid[i] = 0;
    }
    //////////////////////////////////////
    int *cluster_assign = (int *) malloc(sizeof(int) * ndata);
    for (i = 0; i < ndata; i++) {
        cluster_assign[i] = 0;

    }
    /////////////////////////////////////////////////////////////////////////////////////
    /*function call to bisect kmeans*/
    bisecting_kmeans(dim, ndata, data, k, cluster_size, cluster_start, cluster_radius, cluster_centroid,
                     cluster_assign);

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    //function to search the closest point to a randomly generated query point.
    Search(dim, ndata, data, k, cluster_radius, cluster_centroid, cluster_assign, cluster_size, cluster_start);
}