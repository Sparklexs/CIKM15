/************************************************/
/* definitions needed for hash function package */
/************************************************/

#define H1  2197
#define H2  1331
#define HASH(k,i,m) (unsigned int)((( ((long long int) k) * H1 +  \
   ((long long int) i) * (1+(((long long int) k)*H2)%((long long int)(m-1))))\
                          ) % (long long int)(m))
/*
#define HASH(k,i,m) (unsigned int)(( ((long long int) k) * H1 +  \
                       ((long long int) i) * ((1+((long long int) k)*H2)%((long long int)(m-1))))\
                          ) % (long long int)(m)
*/

#include <vector>

using namespace std;


typedef struct {
  int size;
  int * table; 
} hashTable;


// int work = 0; 



/* define here how to get integer key from an element         */
/* NOTE: needs to be adapted to the structure of the elements */
// #define GETKEY(e) (*((int *)(e)))



// hashTable *initHash(int limit, int exact);

hashTable *initHash(int limit, int exact);

void clearHash(hashTable *ht);

void destroyHash(hashTable *ht);

// int insertHash(hashTable *ht, int key, short elem, int overwrite, vector<int> & a, int& hit);

// int insertHash(hashTable *ht, int key, unsigned short elem, int overwrite, vector<int> & a);

// short lookupHash(hashTable *ht, int key, vector<int> & a);

// int GETKEY(unsigned short e, vector<int> & a);




//int insertHash(hashTable *ht, int key, short elem, int overwrite, vector<int> & a);
template<typename T>
int insertHash(hashTable *ht, int key, int elem, int overwrite, vector<T> & a)
{
  // printf("key to insert: %d ", key);
  // printf("elem to insert: %d ", elem);
  int atry;
  unsigned int pos;

  // key = GETKEY(elem, a);
  atry = 0;
  pos = HASH(key, atry, ht->size);
  while ((ht->table[pos] != -1) && (GETKEY(ht->table[pos], a) != key))
  {
    atry++;
    if (atry > 10000)
    {
      printf("Hash-internal error: too many tries during insert!\n");
      return(0);
    }
    pos = HASH(key, atry, ht->size);
  }
    // hit += atry+1;

  if (ht->table[pos] == -1)  
  {
    // printf("pos to insert: %d \n", pos);
    ht->table[pos] = elem;
    return(0);
  }
  else
  { 
    if (overwrite)  ht->table[pos] = elem;
    return pos+1; //The reason +1 is that, pos maybe zero, so all plus one to distinguish from the case above.
    // return(1);
  }
}

template<typename T>
int lookupHash(hashTable *ht, int key, vector<T>& a)
{
  int atry;
  unsigned int pos;

  atry = 0;
  pos = HASH(key, atry, ht->size);
  while ((ht->table[pos] != -1) && (GETKEY(ht->table[pos], a) != key))
  {
    atry++;
    if (atry > 10000)
    {
      printf("Hash-internal error: too many tries during lookup!\n");
      return(0);
    }
    pos = HASH(key, atry, ht->size);
  }

  // printf("The content in the hashtable: %d \n", ht->table[pos]);
  // printf("lookup, find at the pos: %d \n", pos);
  return(ht->table[pos]);
}

template<typename T>
int GETKEY(int e, vector<T> & a){
  if(e<=a.size()){
      // printf("getkey at %d \n", e-1);
      return a[e].did; // store the actual offset
  }else{
    return -1;
  }
}


