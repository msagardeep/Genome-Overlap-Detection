#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

// constants
#define NO_READS 168
#define READ_LENGTH 9000
#define ANCHOR_LENGTH 5
#define NO_ANCHORS READ_LENGTH/ANCHOR_LENGTH
#define BASE 4
#define PRIME 7001

// hash table structure
struct hash
{
	long int val;
	long int pos;
	struct hash *next;
};

// cluster structure
struct cluster
{
	long int read_number;
	long int read_length;
	char sequence[READ_LENGTH];
	long int anchor_vector[NO_ANCHORS];
	long int points[NO_READS];
	long int loc_seq[NO_READS];
	long int loc_cen[NO_READS];
	struct hash *htable[NO_ANCHORS];
	long int counter;
	struct cluster *next;	
};

// initializations
long int total_clusters = 0;
long int yet_to_cluster = NO_READS;
long int cluster_size = 6000;
long int false_collisions = 0;
long int hash_comparisions = 0;
long int base_comparisions = 0;
struct cluster *start = NULL;
struct cluster *last = NULL;

// functions
void pick_cluster_centers(void);
void set_cluster_center(long int, char *);
void cluster_reads(void);
long int determine_cluster(long int, char *);
long int check_subs_matches(struct cluster *, char *, long int, long int, long int);
void print_clusters(void);
void store_clusters(void);

int main(int argc, char **argv)
{
	long int prev_cluster = 0, prev_total_cluster=-1, exit_count = 0;

        printf("Yet to cluster = %ld\n", yet_to_cluster);

        while(yet_to_cluster > 0)
        {
                if(prev_cluster-yet_to_cluster==total_clusters-prev_total_cluster)
                {
                        exit_count++;
                        if(exit_count >= 3*log(NO_READS))
                        {
                                printf("No more clustering ... exiting ...\n\n");
                                store_clusters();
                                return 0;
                        }
                }
                else
                        exit_count = 0;

		hash_comparisions = 0;
		false_collisions = 0;

                prev_cluster = yet_to_cluster;

		prev_total_cluster = total_clusters;

                // copy reads from reads_rem.txt to reads_err.txt
                char sequence[READ_LENGTH] = {0};
                long int read_number;
                FILE *f1, *f2;
                f1 = fopen("reads_rem.txt","r");
                f2 = fopen("reads_err.txt","w");

                while(fscanf(f1,"%ld",&read_number) != EOF)
                {
                        fscanf(f1, "%s", sequence);
                        fprintf(f2, "%ld\n", read_number);
                        fprintf(f2, "%s\n", sequence);
                }
                fclose(f1);
                fclose(f2);

                // pickup random reads and setup the cluster centers
                pick_cluster_centers();

		if(total_clusters > prev_total_cluster)
		{
                	// cluster_reads
                	cluster_reads();

			base_comparisions = hash_comparisions + (false_collisions * 30) + ((prev_cluster-yet_to_cluster) * 30);

                	// print clusters
                	print_clusters();
		}
        }

	store_clusters();

        return 0;
}

void pick_cluster_centers()
{
        long int no_clusters, read_number, count, gap;
        char sequence[READ_LENGTH] = {0};

        // decide number of new clusters
	if(yet_to_cluster > 1)
		no_clusters = log(yet_to_cluster);
	else
		no_clusters = 1;

        count = no_clusters;
	
	last = start;

        // open the given file
        FILE *fp;
        fp = fopen("reads_err.txt","r");

        // set the cluster centers
        while((fscanf(fp,"%ld",&read_number) != EOF) && (count > 0))
        {
                fscanf(fp, "%s", sequence);
               	if(strlen(sequence) >= cluster_size)
               	{
			total_clusters = total_clusters + 1;
                       	set_cluster_center(read_number, sequence);
                       	count = count - 1;
               	}
        }
        fclose(fp);
	cluster_size = 9*cluster_size/10;
}

void set_cluster_center(long int read_number, char *sequence)
{
	struct cluster *new_cluster = NULL;
	long int read_length = strlen(sequence);
	long int no_anchors = read_length/ANCHOR_LENGTH;

	// create a new cluster node and set up the structure fields
	new_cluster = (struct cluster *)malloc(sizeof(struct cluster));

	// set read number
	new_cluster->read_number = read_number;

	// set read_length
	new_cluster->read_length = read_length;

	// store the read
	for(long int i=0; i<read_length; i++)
		new_cluster->sequence[i] = sequence[i];

	// set counter
	new_cluster->counter = -1;

	// set pointer to next cluster as null
	new_cluster->next = NULL;

	// set the hash pointers to null
	for(long int i=0; i<no_anchors; i++)
		new_cluster->htable[i] = NULL;
	
	// compute and store the anchors
	for(long int i=0; i<no_anchors; i++)
	{
		long int anchor_val = 0;
                for(long int j=0; j<ANCHOR_LENGTH; j++)
                {
                        // calculate anchor value
			if(i*ANCHOR_LENGTH + j < read_length)
                       		anchor_val = anchor_val + pow(BASE, j) * (sequence[(i*ANCHOR_LENGTH)+j]-64);
                }
                anchor_val = anchor_val%PRIME;
		new_cluster->anchor_vector[i] = anchor_val;

		// hash the anchor value into htable
		long int hash_val = anchor_val%NO_ANCHORS;
		struct hash *ptr = (struct hash *)malloc(sizeof(struct hash));
		ptr->val = anchor_val;
		ptr->pos = i;
		ptr->next = NULL;
		
		// insert hash node into htable
		ptr->next = new_cluster->htable[hash_val];
		new_cluster->htable[hash_val] = ptr;
	}
	
	// insert new cluster in the list of clusters
	new_cluster->next = start;
	start = new_cluster;
}

void cluster_reads()
{
        long int read_number, remaining = 0;
        char sequence[READ_LENGTH] = {0};

        FILE *f1, *f2;
        f1 = fopen("reads_err.txt", "r");
        f2 = fopen("reads_rem.txt", "w");

        // get the sequences and cluster them
        while(fscanf(f1,"%ld",&read_number) != EOF)
        {
                fscanf(f1, "%s", sequence);
                long int ret = determine_cluster(read_number, sequence);

                // no clustering
                if(ret == 0)
                {
                        remaining = remaining + 1;
                        fprintf(f2, "%ld\n", read_number);
                        fprintf(f2, "%s\n", sequence);
                }
        }

        // these many sequences yet to cluster
        yet_to_cluster = remaining;

        fclose(f1);
        fclose(f2);
}

long int determine_cluster(long int read_number, char *sequence)
{
	long int read_length = strlen(sequence);

	// use rabin karp till 1st actual match is found
	long int anchor_val = 0;
        for(long int i=0; i<ANCHOR_LENGTH; i++)
		if(i<read_length)	
                	anchor_val = anchor_val + pow(BASE,i) * (sequence[i]-64);
        for(long int i=0; i<read_length-ANCHOR_LENGTH+1; i++)
        {
		if(i!=0)
        		anchor_val = (anchor_val - (sequence[i-1]-64)) / 4 
				     + pow(BASE,ANCHOR_LENGTH-1)*(sequence[i+ANCHOR_LENGTH-1]-64);
                anchor_val = anchor_val%PRIME;

		// get the hash value for the anchor
		long int hash_val = anchor_val%NO_ANCHORS;

		// check the anchor value in the clusters
		struct cluster *cluster_ptr = start;

		while(cluster_ptr != last)
		{
			struct hash *hash_ptr = cluster_ptr->htable[hash_val];
			while(hash_ptr != NULL)
			{
				hash_comparisions = hash_comparisions + 1;
				if(hash_ptr->val == anchor_val)
				{
					long int ret_val = 
					   check_subs_matches(cluster_ptr, sequence, i+ANCHOR_LENGTH, hash_ptr->pos,read_length);
					if(ret_val == 1)
					{
						cluster_ptr->counter = cluster_ptr->counter + 1;
						long int counter = cluster_ptr->counter;
						cluster_ptr->points[counter] = read_number;
						cluster_ptr->loc_seq[counter] = i+1;
						cluster_ptr->loc_cen[counter] = hash_ptr->pos*ANCHOR_LENGTH+1;
						return 1;
					}
					false_collisions = false_collisions + 1;
				}
				hash_ptr = hash_ptr->next;
			}
			cluster_ptr = cluster_ptr->next;
		}
        }
	
	return 0;
}

long int check_subs_matches(struct cluster *cluster_ptr, char *sequence, long int begin, long int pos, long int read_length)
{
	for(int i=0; i<30; i++)
	{
		if( (begin+i >= read_length) || (pos+1)*ANCHOR_LENGTH+i >= cluster_ptr->read_length )
			return 0;
		if(sequence[begin+i] != cluster_ptr->sequence[(pos+1)*ANCHOR_LENGTH+i])
			return 0;
	}	
	return 1;
}

void print_clusters()
{
        // print the clusters
	printf("Total reads yet to cluster = %ld\n",yet_to_cluster);
	printf("Total clusters till now = %ld\n",total_clusters);
	printf("Total number of false collisions in this round = %ld\n", false_collisions);
	printf("Total number of comparisions in this round = %ld\n\n", base_comparisions);
}

void store_clusters()
{
	struct cluster *ptr = start;
	FILE *fp;
	char ch1 = '(', ch2 = ')';
	long int count = 0;
	fp = fopen("cluster.txt","a");
	while(ptr != NULL)
	{
		fprintf(fp, "\n%ld\n", ptr->read_number);
		long int read_length = strlen(ptr->sequence);
                for(long int j=0; j<=ptr->counter; j++)
		{
			fprintf(fp, "%c", ch1);
			fprintf(fp, "%ld ", ptr->read_number);
			fprintf(fp, "%ld ", ptr->loc_cen[j]);
			fprintf(fp, "%ld ", ptr->points[j]);
			fprintf(fp, "%ld", ptr->loc_seq[j]);
			fprintf(fp, "%c ", ch2);
		}
		count = count + ptr->counter + 1;
		ptr = ptr->next;
	}
	fclose(fp);
}
