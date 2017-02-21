#include "ruby.h"

static VALUE Kdtree_class;
static VALUE Kdnode_class;

// a node in the tree
typedef struct kdtree_node
{
    float x, y;
    int id;
    struct kdtree_node *left;
    struct kdtree_node *right;
} kdtree_node;

// the tree itself
typedef struct kdtree
{
    struct kdtree_node *root;
} kdtree;

// helper macro for digging out our structs
#define KDTREEP \
    struct kdtree *kdtreep; \
    Data_Get_Struct(kdtree, struct kdtree, kdtreep);

#define KDNODEP \
    struct kdtree_node *node; \
    Data_Get_Struct(kdnode, struct kdtree_node, node);

static VALUE kdtree_node_alloc(VALUE klass);
static void kdtree_node_free(struct kdtree_node *node);
static VALUE kdtree_node_initialize(VALUE kdnode, VALUE x, VALUE y, VALUE data);
static VALUE kdtree_node_id(VALUE kdnode);
static VALUE kdtree_node_x(VALUE kdnode);
static VALUE kdtree_node_y(VALUE kdnode);

static VALUE kdtree_alloc(VALUE klass);
static void kdtree_free(struct kdtree *kdtreep);
static void kdtree_free_node(struct kdtree_node *node);
static VALUE kdtree_initialize(VALUE kdtree, VALUE points);
static kdtree_node* kdtree_build(struct kdtree_node *nodes, int min, int max, int depth);

// Heap implementation

#define HEAP_LCHILD(x) 2 * x + 1
#define HEAP_RCHILD(x) 2 * x + 2
#define HEAP_PARENT(x) (x - 1) / 2

typedef struct heapNode {
    float distance;
    kdtree_node *node;
} heapNode;

typedef struct maxHeap {
    int len;
    int size;
    heapNode *elem;
} maxHeap;


static maxHeap initmaxHeap() {
    maxHeap hp ;
    hp.size = 0 ;
    return hp ;
}

static void deletemaxHeap(maxHeap *hp) {
    free(hp->elem) ;
}

static void swap(heapNode *n1, heapNode *n2) {
    heapNode temp = *n1 ;
    *n1 = *n2 ;
    *n2 = temp ;
}

static void heapify(maxHeap *hp, int i) {
  int largest = (HEAP_LCHILD(i) < hp->len && hp->elem[HEAP_LCHILD(i)].distance > hp->elem[i].distance) ? HEAP_LCHILD(i) : i ;
  if(HEAP_RCHILD(i) < hp->len && hp->elem[HEAP_RCHILD(i)].distance > hp->elem[largest].distance) {
    largest = HEAP_RCHILD(i) ;
  }
  if(largest != i) {
    swap(&(hp->elem[i]), &(hp->elem[largest])) ;
    heapify(hp, largest) ;
  }
}

static kdtree_node *pop(maxHeap *hp) {
  if (!hp->len) {
    return NULL;
  }
  kdtree_node *node = hp->elem[0].node;
  hp->elem[0] = hp->elem[--(hp->len)];
  heapify(hp, 0);
  return node;
}

static void push(maxHeap *hp, float distance, kdtree_node *node, int maxNodes) {
  if (hp->len + 1 >= hp->size) {
    hp->size = hp->size ? hp->size * 2 : 4;
    hp->elem = realloc(hp->elem, hp->size * sizeof(heapNode));
  }

  heapNode nd;
  nd.node = node;
  nd.distance = distance;

  int i = (hp->len)++;
  while(i && nd.distance > hp->elem[HEAP_PARENT(i)].distance) {
        hp->elem[i] = hp->elem[HEAP_PARENT(i)] ;
        i = HEAP_PARENT(i) ;
    }
  // while(i && nd.distance < hp->elem[HEAP_PARENT(i)].distance) {
  //     hp->elem[i] = hp->elem[HEAP_PARENT(i)] ;
  //     i = HEAP_PARENT(i);
  // }
  hp->elem[i] = nd;
  if (hp->len > maxNodes) {
    pop(hp);
  }
}

static heapNode *peek(maxHeap *hp) {
  return &hp->elem[0];
}

// Node

static VALUE kdtree_node_alloc(VALUE klass)
{
  struct kdtree_node *node;
  VALUE obj = Data_Make_Struct(klass, struct kdtree_node, 0, kdtree_node_free, node);
  return obj;
}

static VALUE kdtree_node_initialize(VALUE kdnode, VALUE x, VALUE y, VALUE id)
{
  KDNODEP;
  node->id = NUM2INT(id);
  node->x = NUM2DBL(x);
  node->y = NUM2DBL(y);
  return kdnode;
}

static void kdtree_node_free(struct kdtree_node *node)
{
}

static VALUE kdtree_node_id(VALUE kdnode)
{
  KDNODEP;
  return INT2NUM(node->id);
}

static VALUE kdtree_node_x(VALUE kdnode)
{
  KDNODEP;
  return DBL2NUM(node->x);
}

static VALUE kdtree_node_y(VALUE kdnode)
{
  KDNODEP;
  return DBL2NUM(node->y);
}

static VALUE kdtree_node_left(VALUE kdnode)
{
  KDNODEP;

  if (!node->left) {
    return Qnil;
  }

  return Data_Wrap_Struct(Kdnode_class, NULL, kdtree_node_free, node->left);
}

static VALUE kdtree_node_right(VALUE kdnode)
{
  KDNODEP;

  if (!node->right) {
    return Qnil;
  }

  return Data_Wrap_Struct(Kdnode_class, NULL, kdtree_node_free, node->right);
}

static VALUE kdtree_node_to_s(VALUE kdnode)
{
    char buf[256];
    KDNODEP;

    sprintf(buf, "#<%s:%p id=%d (%f,%f)>", rb_obj_classname(kdnode), (void*)kdnode, node->id, node->x, node->y);
    return rb_str_new(buf, strlen(buf));
}

// Tree

static VALUE kdtree_alloc(VALUE klass)
{
    struct kdtree *kdtreep;
    VALUE obj = Data_Make_Struct(klass, struct kdtree, 0, kdtree_free, kdtreep);
    kdtreep->root = NULL;
    return obj;
}

static void kdtree_free(struct kdtree *kdtreep)
{
    if (kdtreep && kdtreep->root) {
        kdtree_free_node(kdtreep->root);
    }
}

static void kdtree_free_node(struct kdtree_node *node)
{
  // if (node) {
  //   if (node->left) {
  //     kdtree_free_node(node->left);
  //   }
  //   if (node->right) {
  //     kdtree_free_node(node->right);
  //   }
  //   free(node);
  // }
}

static VALUE kdtree_initialize(VALUE kdtree, VALUE points) {
  KDTREEP;

  if (TYPE(points) == T_ARRAY) {
    int i;
    int len = RARRAY_LEN(points);
    struct kdtree_node *nodes = ALLOC_N(struct kdtree_node, len);
    for (i = 0; i < RARRAY_LEN(points); ++i) {
      VALUE ptr = rb_ary_entry(points, i);
      struct kdtree_node *node;
      Data_Get_Struct(ptr, struct kdtree_node, node);
      nodes[i] = *node;
    }
    kdtreep->root = kdtree_build(nodes, 0, len, 0);
  } else {
    rb_raise(rb_eTypeError, "array required to init KDTree");
  }
  return kdtree;
}

static VALUE kdtree_root(VALUE kdtree) {
  KDTREEP;

  if (!kdtreep->root) {
    return Qnil;
  }

  return Data_Wrap_Struct(Kdnode_class, NULL, kdtree_node_free, kdtreep->root);
}

static int comparex(const void *pa, const void *pb)
{
    float a = ((const struct kdtree_node*)pa)->x;
    float b = ((const struct kdtree_node*)pb)->x;
    return (a < b) ? -1 : ((a > b) ? 1 : 0);
}

static int comparey(const void *pa, const void *pb)
{
    float a = ((const struct kdtree_node*)pa)->y;
    float b = ((const struct kdtree_node*)pb)->y;
    return (a < b) ? -1 : ((a > b) ? 1 : 0);
}

static float kdnode_get_dim(struct kdtree_node *node, int depth)
{
  if (depth % 2) {
    return node->y;
  } else {
    return node->x;
  }
}

static kdtree_node* kdtree_build(struct kdtree_node *nodes, int min, int max, int depth)
{
  if (max <= min) {
    return NULL;
  }

  struct kdtree_node *m;
  int median;

  int(*compare)(const void *, const void *);
  compare = (depth % 2) ? comparex : comparey;

  qsort(nodes + min, max - min, sizeof(struct kdtree_node), compare);

  median = (min + max) / 2;
  m = nodes + median;
  m->left = kdtree_build(nodes, min, median, depth + 1);
  m->right = kdtree_build(nodes, median + 1, max, depth + 1);
  return m;
}

static void nearest_search(maxHeap *best, float maxNodes, float x, float y, kdtree_node *node, int depth)
{
  if (!node) {
    return;
  }

  float dist = pow(x - node->x, 2) + pow(y - node->y, 2);
  float diff;

  if (depth % 2) {
    diff = x - node->x;
  } else {
    diff = y - node->y;
  }

  if (best->len < maxNodes || dist <= peek(best)->distance) {
    push(best, dist, node, maxNodes);
  }

  kdtree_node *near = diff <= 0 ? node->left : node->right;
  kdtree_node *far = diff <= 0 ? node->right : node->left;

  nearest_search(best, maxNodes, x, y, near, depth+1);

  if (best->len < maxNodes || (diff * diff) < peek(best)->distance) {
    nearest_search(best, maxNodes, x, y, far, depth+1);
  }
}

static kdtree_node* kdtree_insert(struct kdtree_node *p, struct kdtree_node *n, int depth)
{
    int(*compare)(const void *, const void *);
    compare = (depth % 2) ? comparex : comparey;
    if (n == NULL) {
        n = p;
    } else if (n == p) {
        // nothing
    } else {
        if (compare(p, n) == -1) {
            n->left = kdtree_insert(p, n->left, depth + 1);
        } else {
            n->right = kdtree_insert(p, n->right, depth + 1);
        }
    }
    return n;
}

static kdtree_node* kdtree_find_min(struct kdtree_node *node, int dim, int depth)
{
  if (node == NULL) {
    return NULL;
  }

  if (depth %2 == dim) {
    if (node->left != NULL) {
      return kdtree_find_min(node->left, dim, depth+1);
    }
    return node;
  }

  struct kdtree_node *min = node;

  struct kdtree_node *left = kdtree_find_min(node->left, dim, depth+1);
  struct kdtree_node *right = kdtree_find_min(node->right, dim, depth+1);

  int(*compare)(const void *, const void *);
  compare = (dim % 2) ? comparex : comparey;

  if (left != NULL && compare(left, min) == -1) {
    min = left;
  }
  if (right != NULL && compare(right, min) == -1) {
    min = right;
  }

  return min;
}

static void kdtree_copy(struct kdtree_node *dest, struct kdtree_node *source) {
  dest->id = source->id;
  dest->x = source->x;
  dest->y = source->y;
}

static kdtree_node* kdtree_remove(struct kdtree_node *p, struct kdtree_node *node, int depth)
{
  if (p == NULL) {
    return NULL;
  }

  int dim = depth % 2;

  if (p->id == node->id) {
    if (p->right != NULL) {
      struct kdtree_node *min = kdtree_find_min(p->right, dim, depth + 1);
      kdtree_copy(p, min);
      p->right = kdtree_remove(p->right, min, depth + 1);
    } else if (p->left != NULL) {
      struct kdtree_node *min = kdtree_find_min(p->left, dim, depth + 1);
      kdtree_copy(p, min);
      p->right = kdtree_remove(p->left, min, depth + 1);
      p->left = NULL;
    } else {
      //free(p);
      return NULL;
    }
    return p;
  }

  int(*compar)(const void *, const void *);
  compar = (dim % 2) ? comparex : comparey;

  if (compar(node, p) == -1) {
    p->left = kdtree_remove(p->left, node, depth + 1);
  } else {
    p->right = kdtree_remove(p->right, node, depth + 1);
  }

  return p;
}

static VALUE kdtree_add(VALUE kdtree, VALUE kdnode) {
    KDTREEP;
    KDNODEP;

    kdtreep->root = kdtree_insert(node, kdtreep->root, 0);
    return kdtree;
}

static VALUE kdtree_delete(VALUE kdtree, VALUE kdnode) {
  KDTREEP;
  KDNODEP;

  kdtreep->root = kdtree_remove(kdtreep->root, node, 0);

  return kdtree;
}

static VALUE kdtree_nearest(int argc, VALUE* argv, VALUE kdtree)
{
  KDTREEP;
  float x = NUM2DBL(argv[0]);
  float y = NUM2DBL(argv[1]);
  int maxNodes = 1;

  if (argc > 2) {
    maxNodes = NUM2INT(argv[2]);
  }

  maxHeap *h = (maxHeap *)calloc(1, sizeof (maxHeap));

  if (argc > 3) {
    float maxDistance = NUM2DBL(argv[3]);
    for (int i=0; i<maxNodes; ++i) {
      push(h, maxDistance, NULL, maxNodes);
    }
  }

  if (kdtreep->root) {
    nearest_search(h, maxNodes, x, y, kdtreep->root, 0);
  }


  if (argc <= 2) {
    // no maxNodes, return one result
    return Data_Wrap_Struct(Kdnode_class, NULL, kdtree_node_free, pop(h));
  } else {
    // return array
    VALUE result = rb_ary_new();
    int cnt = h->len;
    if (maxNodes < cnt) {
      cnt = maxNodes;
    }

    struct kdtree_node *node;

    for (int i=0; i<cnt; ++i) {
      node = pop(h);
      if (node) {
        VALUE v = Data_Wrap_Struct(Kdnode_class, NULL, kdtree_node_free, node);
        rb_ary_push(result, v);
      }
    }
    return result;
  }
}

void Init_kdtree()
{
    Kdnode_class = rb_define_class("KDTreeNode", rb_cObject);
    Kdtree_class = rb_define_class("KDTree", rb_cObject);

    rb_define_alloc_func(Kdtree_class, kdtree_alloc);
    rb_define_method(Kdtree_class, "initialize", kdtree_initialize, 1);
    rb_define_method(Kdtree_class, "nearest", kdtree_nearest, -1);
    rb_define_method(Kdtree_class, "<<", kdtree_add, 1);
    rb_define_method(Kdtree_class, "delete", kdtree_delete, 1);
    rb_define_method(Kdtree_class, "root", kdtree_root, 0);

    rb_define_alloc_func(Kdnode_class, kdtree_node_alloc);
    rb_define_method(Kdnode_class, "initialize", kdtree_node_initialize, 3);
    rb_define_method(Kdnode_class, "id", kdtree_node_id, 0);
    rb_define_method(Kdnode_class, "x", kdtree_node_x, 0);
    rb_define_method(Kdnode_class, "y", kdtree_node_y, 0);
    rb_define_method(Kdnode_class, "to_s", kdtree_node_to_s, 0);
    rb_define_method(Kdnode_class, "inspect", kdtree_node_to_s, 0);
    rb_define_method(Kdnode_class, "left", kdtree_node_left, 0);
    rb_define_method(Kdnode_class, "right", kdtree_node_right, 0);
}
