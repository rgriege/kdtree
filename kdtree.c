/*
This file is part of ``kdtree'', a library for working with kd-trees.
Copyright (C) 2007-2011 John Tsiombikas <nuclear@member.fsf.org>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.
*/
/* single nearest neighbor search written by Tamas Nepusz <tamas@cs.rhul.ac.uk> */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "kdtree.h"

#if defined(WIN32) || defined(__WIN32__)
#include <malloc.h>
#endif

#ifdef USE_LIST_NODE_ALLOCATOR

#ifndef NO_PTHREADS
#include <pthread.h>
#else

#ifndef I_WANT_THREAD_BUGS
#error "You are compiling with the fast list node allocator, with pthreads disabled! This WILL break if used from multiple threads."
#endif	/* I want thread bugs */

#endif	/* pthread support */
#endif	/* use list node allocator */

struct kdhyperrect {
	int dim;
	kdcoord minmax[];
};

struct kdnode {
	struct kdnode *left, *right;	/* negative/positive side */
	void *data;
	int dir;
	kdcoord pos[];
};

struct res_node {
	struct kdnode *item;
	kdcoord dist_sq;
	struct res_node *next;
};

struct kdtree {
	int dim;
	int dim_zero_filter;
	struct kdnode *root;
	struct kdhyperrect *rect;
	struct kdhyperrect *rect_copy;
	void (*destr)(void*);
};

struct kdres {
	struct kdtree *tree;
	struct res_node *rlist, *riter;
	int size;
};

#define SQ(x)			((x) * (x))

static struct kdnode* node_create(struct kdtree *tree, const kdcoord *pos, void *data);

static void clear_rec(struct kdnode *node, void (*destr)(void*));
static void insert_rec(struct kdnode **nptr, struct kdnode *node, int dir, int dim);
static int rlist_insert(struct res_node *list, struct kdnode *item, kdcoord dist_sq);
static void clear_results(struct kdres *set);

static size_t hyperrect_size(int dim);
static struct kdhyperrect* hyperrect_create_in_buffer(int dim, void *mem);
static void hyperrect_clear(struct kdhyperrect *rect);
static kdcoord *hyperrect_min(struct kdhyperrect *rect);
static kdcoord *hyperrect_max(struct kdhyperrect *rect);
static void hyperrect_copy(struct kdhyperrect *rect, const kdcoord *min, const kdcoord *max);
static void hyperrect_extend(struct kdhyperrect *rect, const kdcoord *pos);
static kdcoord hyperrect_dist_sq(struct kdhyperrect *rect, const kdcoord *pos);

#ifdef USE_LIST_NODE_ALLOCATOR
static struct res_node *alloc_resnode(void);
static void free_resnode(struct res_node*);
#else
#define alloc_resnode()		malloc(sizeof(struct res_node))
#define free_resnode(n)		free(n)
#endif


size_t kd_tree_size(int k)
{
	size_t tree_size = sizeof(struct kdtree);
	size_t rect_size = hyperrect_size(k);
	return tree_size + 2 * rect_size;
}

size_t kd_node_size(struct kdtree *tree)
{
	return sizeof(struct kdnode) + tree->dim * sizeof(kdcoord);
}

struct kdtree *kd_create(int k, int k_zero_filter)
{
	size_t size = kd_tree_size(k);
	void *mem;

	if(!(mem = malloc(size))) {
		return 0;
	}

	return kd_create_in_buffer(k, k_zero_filter, mem);
}

struct kdtree *kd_create_in_buffer(int k, int k_zero_filter, void *mem)
{
	size_t tree_size = sizeof(struct kdtree);
	size_t rect_size = hyperrect_size(k);
	struct kdtree *tree = mem;

	tree->dim = k;
	tree->dim_zero_filter = k_zero_filter;
	tree->root = 0;
	tree->destr = 0;
	tree->rect = hyperrect_create_in_buffer(tree->dim, (char*)mem + tree_size);
	tree->rect_copy = hyperrect_create_in_buffer(tree->dim, (char*)mem + tree_size + rect_size);

	return tree;
}

void kd_free(struct kdtree *tree)
{
	if(tree) {
		kd_clear(tree);
		tree->rect = 0;
		tree->rect_copy = 0;
		free(tree);
	}
}

static void clear_rec(struct kdnode *node, void (*destr)(void*))
{
	if(!node) return;

	clear_rec(node->left, destr);
	clear_rec(node->right, destr);
	
	if(destr) {
		destr(node->data);
	}
	free(node);
}

void kd_clear(struct kdtree *tree)
{
	clear_rec(tree->root, tree->destr);
	tree->root = 0;

	hyperrect_clear(tree->rect);
	hyperrect_clear(tree->rect_copy);
}

void kd_data_destructor(struct kdtree *tree, void (*destr)(void*))
{
	tree->destr = destr;
}

static struct kdnode* node_create(struct kdtree *tree, const kdcoord *pos, void *data)
{
	size_t size = kd_node_size(tree);
	void *mem;

	if(!(mem = malloc(size))) {
		return 0;
	}
	return kd_create_node_in_buffer(tree, pos, data, mem);
}

static void insert_rec(struct kdnode **nptr, struct kdnode *node, int dir, int dim)
{
	int new_dir;
	struct kdnode *branch;

	if(!*nptr) {
		node->dir = dir;
		*nptr = node;
		return;
	}

	branch = *nptr;
	new_dir = (branch->dir + 1) % dim;
	if(node->pos[branch->dir] < branch->pos[branch->dir]) {
		insert_rec(&branch->left, node, new_dir, dim);
	} else {
		insert_rec(&branch->right, node, new_dir, dim);
	}
}

struct kdnode *kd_create_node_in_buffer(struct kdtree *tree, const kdcoord *pos, void *data, void *mem)
{
	struct kdnode *node = mem;

	memcpy(node->pos, pos, tree->dim * sizeof *pos);
	node->data = data;
	node->dir = 0;
	node->left = node->right = 0;

	return node;
}

void kd_insert_node(struct kdtree *tree, struct kdnode *node)
{
	int was_empty = tree->root == 0;

	insert_rec(&tree->root, node, 0, tree->dim);

	if (was_empty) {
		hyperrect_copy(tree->rect, node->pos, node->pos);
	} else {
		hyperrect_extend(tree->rect, node->pos);
	}
}

int kd_insert(struct kdtree *tree, const kdcoord *pos, void *data)
{
	struct kdnode *node;

	if(!(node = node_create(tree, pos, data))) {
		return -1;
	}

	kd_insert_node(tree, node);
	return 0;
}

int kd_insert3(struct kdtree *tree, kdcoord x, kdcoord y, kdcoord z, void *data)
{
	kdcoord buf[3];
	buf[0] = x;
	buf[1] = y;
	buf[2] = z;
	return kd_insert(tree, buf, data);
}

struct kdnode *kd_root(struct kdtree *tree)
{
	return tree->root;
}

struct kdnode *kd_node_left(struct kdnode *node)
{
	return node->left;
}

struct kdnode *kd_node_right(struct kdnode *node)
{
	return node->right;
}

static int find_nearest(struct kdnode *node, const kdcoord *pos, kdcoord range, struct res_node *list, int ordered, int dim)
{
	kdcoord dist_sq, dx;
	int i, ret, added_res = 0;

	if(!node) return 0;

	dist_sq = 0;
	for(i=0; i<dim; i++) {
		dist_sq += SQ(node->pos[i] - pos[i]);
	}
	if(dist_sq <= SQ(range)) {
		if(rlist_insert(list, node, ordered ? dist_sq : -1.0) == -1) {
			return -1;
		}
		added_res = 1;
	}

	dx = pos[node->dir] - node->pos[node->dir];

	ret = find_nearest(dx <= 0.0 ? node->left : node->right, pos, range, list, ordered, dim);
	if(ret >= 0 && fabs(dx) < range) {
		added_res += ret;
		ret = find_nearest(dx <= 0.0 ? node->right : node->left, pos, range, list, ordered, dim);
	}
	if(ret == -1) {
		return -1;
	}
	added_res += ret;

	return added_res;
}

static int kd_filter_zero_weight(const kdcoord *goal, const kdcoord *node, int dim)
{
	int i;
	for(i=0; i < dim; i++)
		if (goal[i] == 0.0 && node[i] > 0.0)
			return -1;
	return 0;
}

static void kd_nearest_i(struct kdnode *node, const kdcoord *pos, int dim_zero_filter, struct kdnode **result, kdcoord *result_dist_sq, struct kdhyperrect* rect)
{
	int dir = node->dir;
	int i;
	kdcoord dummy, dist_sq;
	struct kdnode *nearer_subtree, *farther_subtree;
	kdcoord *nearer_hyperrect_coord, *farther_hyperrect_coord;

	/* Decide whether to go left or right in the tree */
	dummy = pos[dir] - node->pos[dir];
	if (dummy <= 0) {
		nearer_subtree = node->left;
		farther_subtree = node->right;
		nearer_hyperrect_coord = hyperrect_max(rect) + dir;
		farther_hyperrect_coord = hyperrect_min(rect) + dir;
	} else {
		nearer_subtree = node->right;
		farther_subtree = node->left;
		nearer_hyperrect_coord = hyperrect_min(rect) + dir;
		farther_hyperrect_coord = hyperrect_max(rect) + dir;
	}

	if (nearer_subtree) {
		/* Slice the hyperrect to get the hyperrect of the nearer subtree */
		dummy = *nearer_hyperrect_coord;
		*nearer_hyperrect_coord = node->pos[dir];
		/* Recurse down into nearer subtree */
		kd_nearest_i(nearer_subtree, pos, dim_zero_filter, result, result_dist_sq, rect);
		/* Undo the slice */
		*nearer_hyperrect_coord = dummy;
	}

	/* Check the distance of the point at the current node, compare it
	 * with our best so far */
	if (kd_filter_zero_weight(pos, node->pos, dim_zero_filter) == 0) {
		dist_sq = 0;
		for(i=0; i < rect->dim; i++)
			dist_sq += SQ(node->pos[i] - pos[i]);
		if (dist_sq < *result_dist_sq) {
			*result = node;
			*result_dist_sq = dist_sq;
		}
	}

	if (farther_subtree) {
		/* Get the hyperrect of the farther subtree */
		dummy = *farther_hyperrect_coord;
		*farther_hyperrect_coord = node->pos[dir];
		/* Check if we have to recurse down by calculating the closest
		 * point of the hyperrect and see if it's closer than our
		 * minimum distance in result_dist_sq. */
		if (hyperrect_dist_sq(rect, pos) < *result_dist_sq) {
			/* Recurse down into farther subtree */
			kd_nearest_i(farther_subtree, pos, dim_zero_filter, result, result_dist_sq, rect);
		}
		/* Undo the slice on the hyperrect */
		*farther_hyperrect_coord = dummy;
	}
}

static struct kdnode *kd_nearest_1(struct kdtree *kd, const kdcoord *pos)
{
	struct kdhyperrect *rect;
	struct kdnode *result;
	kdcoord dist_sq;
	int i;

	if (!kd) return 0;
	if (!kd->root) return 0;

	/* Duplicate the bounding hyperrectangle, we will work on the copy.
	 * This copy is kept around permanently to avoid dynamic allocations. */
	rect = kd->rect_copy;
	hyperrect_copy(rect, hyperrect_min(kd->rect), hyperrect_max(kd->rect));

	/* Our first guesstimate is the root node */
	if (kd_filter_zero_weight(pos, kd->root->pos, kd->dim_zero_filter) == 0) {
		result = kd->root;
		dist_sq = 0;
		for (i = 0; i < kd->dim; i++)
			dist_sq += SQ(result->pos[i] - pos[i]);
	} else {
		result = NULL;
#ifdef KDTREE_USE_FLOAT
		dist_sq = FLT_MAX;
#else
		dist_sq = DBL_MAX;
#endif
	}

	/* Search for the nearest neighbour recursively */
	kd_nearest_i(kd->root, pos, kd->dim_zero_filter, &result, &dist_sq, rect);

	return result;
}

struct kdres *kd_nearest(struct kdtree *kd, const kdcoord *pos)
{
	struct kdnode *result;
	struct kdres *rset;

	/* Allocate result set */
	if(!(rset = malloc(sizeof *rset))) {
		return 0;
	}
	if(!(rset->rlist = alloc_resnode())) {
		free(rset);
		return 0;
	}
	rset->rlist->next = 0;
	rset->tree = kd;

	result = kd_nearest_1(kd, pos);

	/* Store the result */
	if (result) {
		if (rlist_insert(rset->rlist, result, -1.0) == -1) {
			kd_res_free(rset);
			return 0;
		}
		rset->size = 1;
		kd_res_rewind(rset);
		return rset;
	} else {
		kd_res_free(rset);
		return 0;
	}
}


struct kdres *kd_nearest3(struct kdtree *tree, kdcoord x, kdcoord y, kdcoord z)
{
	kdcoord pos[3];
	pos[0] = x;
	pos[1] = y;
	pos[2] = z;
	return kd_nearest(tree, pos);
}

int kd_nearest_one(struct kdtree *kd, const kdcoord *pos, kdcoord *out, void **data)
{
	struct kdnode *result = kd_nearest_1(kd, pos);
	if(!result) {
		return 0;
	}
	if(out) {
		memcpy(out, result->pos, kd->dim * sizeof *pos);
	}
	if(data) {
		*data = result->data;
	}
	return 1;
}

struct kdres *kd_nearest_range(struct kdtree *kd, const kdcoord *pos, kdcoord range)
{
	int ret;
	struct kdres *rset;

	if(!(rset = malloc(sizeof *rset))) {
		return 0;
	}
	if(!(rset->rlist = alloc_resnode())) {
		free(rset);
		return 0;
	}
	rset->rlist->next = 0;
	rset->tree = kd;

	if((ret = find_nearest(kd->root, pos, range, rset->rlist, 0, kd->dim)) == -1) {
		kd_res_free(rset);
		return 0;
	}
	rset->size = ret;
	kd_res_rewind(rset);
	return rset;
}

struct kdres *kd_nearest_range3(struct kdtree *tree, kdcoord x, kdcoord y, kdcoord z, kdcoord range)
{
	kdcoord buf[3];
	buf[0] = x;
	buf[1] = y;
	buf[2] = z;
	return kd_nearest_range(tree, buf, range);
}

void kd_res_free(struct kdres *rset)
{
	clear_results(rset);
	free_resnode(rset->rlist);
	free(rset);
}

int kd_res_size(struct kdres *set)
{
	return (set->size);
}

void kd_res_rewind(struct kdres *rset)
{
	rset->riter = rset->rlist->next;
}

int kd_res_end(struct kdres *rset)
{
	return rset->riter == 0;
}

int kd_res_next(struct kdres *rset)
{
	rset->riter = rset->riter->next;
	return rset->riter != 0;
}

void *kd_res_item(struct kdres *rset, kdcoord *pos)
{
	if(rset->riter) {
		if(pos) {
			memcpy(pos, rset->riter->item->pos, rset->tree->dim * sizeof *pos);
		}
		return rset->riter->item->data;
	}
	return 0;
}

void *kd_res_item3(struct kdres *rset, kdcoord *x, kdcoord *y, kdcoord *z)
{
	if(rset->riter) {
		if(x) *x = rset->riter->item->pos[0];
		if(y) *y = rset->riter->item->pos[1];
		if(z) *z = rset->riter->item->pos[2];
		return rset->riter->item->data;
	}
	return 0;
}

void *kd_res_item_data(struct kdres *set)
{
	return kd_res_item(set, 0);
}

/* ---- hyperrectangle helpers ---- */
static size_t hyperrect_size(int dim)
{
	return sizeof(struct kdhyperrect) + 2 * dim * sizeof(kdcoord);
}

static struct kdhyperrect* hyperrect_create_in_buffer(int dim, void *mem)
{
	struct kdhyperrect* rect = mem;

	rect->dim = dim;
	hyperrect_clear(rect);

	return rect;
}

static void hyperrect_clear(struct kdhyperrect *rect)
{
	size_t size = 2 * rect->dim * sizeof(kdcoord);
	memset(rect->minmax, 0, size);
}

static void hyperrect_copy(struct kdhyperrect *rect, const kdcoord *min, const kdcoord *max)
{
	size_t size = rect->dim * sizeof(kdcoord);
	kdcoord *pmin = hyperrect_min(rect);
	kdcoord *pmax = hyperrect_max(rect);
	memcpy(pmin, min, size);
	memcpy(pmax, max, size);
}

static kdcoord *hyperrect_min(struct kdhyperrect *rect)
{
	return rect->minmax;
}

static kdcoord *hyperrect_max(struct kdhyperrect *rect)
{
	return rect->minmax + rect->dim;
}

static void hyperrect_extend(struct kdhyperrect *rect, const kdcoord *pos)
{
	kdcoord *min = hyperrect_min(rect);
	kdcoord *max = hyperrect_max(rect);
	int i;

	for (i=0; i < rect->dim; i++) {
		if (pos[i] < min[i]) {
			min[i] = pos[i];
		}
		if (pos[i] > max[i]) {
			max[i] = pos[i];
		}
	}
}

static kdcoord hyperrect_dist_sq(struct kdhyperrect *rect, const kdcoord *pos)
{
	kdcoord *min = hyperrect_min(rect);
	kdcoord *max = hyperrect_max(rect);
	int i;
	kdcoord result = 0;

	for (i=0; i < rect->dim; i++) {
		if (pos[i] < min[i]) {
			result += SQ(min[i] - pos[i]);
		} else if (pos[i] > max[i]) {
			result += SQ(max[i] - pos[i]);
		}
	}

	return result;
}

/* ---- static helpers ---- */

#ifdef USE_LIST_NODE_ALLOCATOR
/* special list node allocators. */
static struct res_node *free_nodes;

#ifndef NO_PTHREADS
static pthread_mutex_t alloc_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

static struct res_node *alloc_resnode(void)
{
	struct res_node *node;

#ifndef NO_PTHREADS
	pthread_mutex_lock(&alloc_mutex);
#endif

	if(!free_nodes) {
		node = malloc(sizeof *node);
	} else {
		node = free_nodes;
		free_nodes = free_nodes->next;
		node->next = 0;
	}

#ifndef NO_PTHREADS
	pthread_mutex_unlock(&alloc_mutex);
#endif

	return node;
}

static void free_resnode(struct res_node *node)
{
#ifndef NO_PTHREADS
	pthread_mutex_lock(&alloc_mutex);
#endif

	node->next = free_nodes;
	free_nodes = node;

#ifndef NO_PTHREADS
	pthread_mutex_unlock(&alloc_mutex);
#endif
}
#endif	/* list node allocator or not */


/* inserts the item. if dist_sq is >= 0, then do an ordered insert */
/* TODO make the ordering code use heapsort */
static int rlist_insert(struct res_node *list, struct kdnode *item, kdcoord dist_sq)
{
	struct res_node *rnode;

	if(!(rnode = alloc_resnode())) {
		return -1;
	}
	rnode->item = item;
	rnode->dist_sq = dist_sq;

	if(dist_sq >= 0.0) {
		while(list->next && list->next->dist_sq < dist_sq) {
			list = list->next;
		}
	}
	rnode->next = list->next;
	list->next = rnode;
	return 0;
}

static void clear_results(struct kdres *rset)
{
	struct res_node *tmp, *node = rset->rlist->next;

	while(node) {
		tmp = node;
		node = node->next;
		free_resnode(tmp);
	}

	rset->rlist->next = 0;
}
