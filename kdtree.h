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
#ifndef _KDTREE_H_
#define _KDTREE_H_

#ifdef KDTREE_USE_FLOAT
#define kdcoord float
#else
#define kdcoord double
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct kdtree;
struct kdres;
struct kdnode;

/* compute memory buffer sizes for the tree and nodes.
 * if you use these to manually allocated the tree and nodes,
 * do NOT call kd_free or kd_clear */
size_t kd_tree_size(int k);
size_t kd_node_size(struct kdtree *tree);

/* create a kd-tree for "k"-dimensional data.
 * if the k_zero_filter parameter is not 0, then any nearest query
 * that has a coord of 0 in any position will filter out any nodes
 * with non-0 values in those positions */
struct kdtree *kd_create(int k, int k_zero_filter);
struct kdtree *kd_create_in_buffer(int k, int k_zero_filter, void *mem);

/* free the struct kdtree */
void kd_free(struct kdtree *tree);

/* remove all the elements from the tree */
void kd_clear(struct kdtree *tree);

/* if called with non-null 2nd argument, the function provided
 * will be called on data pointers (see kd_insert) when nodes
 * are to be removed from the tree.
 */
void kd_data_destructor(struct kdtree *tree, void (*destr)(void*));

/* insert a node, specifying its position, and optional data */
int kd_insert(struct kdtree *tree, const kdcoord *pos, void *data);
int kd_insert3(struct kdtree *tree, kdcoord x, kdcoord y, kdcoord z, void *data);

struct kdnode *kd_create_node_in_buffer(struct kdtree *tree, const kdcoord *pos, void *data, void *mem);
void kd_insert_node(struct kdtree *tree, struct kdnode *node);

/* node navigation */
struct kdnode *kd_root(struct kdtree *tree);
struct kdnode *kd_node_left(struct kdnode *node);
struct kdnode *kd_node_right(struct kdnode *node);

/* Find the nearest node from a given point.
 *
 * This function returns a pointer to a result set with at most one element.
 */
struct kdres *kd_nearest(struct kdtree *tree, const kdcoord *pos);
struct kdres *kd_nearest3(struct kdtree *tree, kdcoord x, kdcoord y, kdcoord z);

/* Find the nearest node from a given point.
 *
 * This function returns non-zero on success or zero if a node was
 * not found. On success, the pointer to the data of the nearest node
 * is also set. Unlike kd_nearest*, it does not perform any allocations
 * for the result set.
 */
int kd_nearest_one(struct kdtree *tree, const kdcoord *pos, kdcoord *out, void **data);

/* Find any nearest nodes from a given point within a range.
 *
 * This function returns a pointer to a result set, which can be manipulated
 * by the kd_res_* functions.
 * The returned pointer can be null as an indication of an error. Otherwise
 * a valid result set is always returned which may contain 0 or more elements.
 * The result set must be deallocated with kd_res_free after use.
 */
struct kdres *kd_nearest_range(struct kdtree *tree, const kdcoord *pos, kdcoord range);
struct kdres *kd_nearest_range3(struct kdtree *tree, kdcoord x, kdcoord y, kdcoord z, kdcoord range);

/* frees a result set returned by kd_nearest_range() */
void kd_res_free(struct kdres *set);

/* returns the size of the result set (in elements) */
int kd_res_size(struct kdres *set);

/* rewinds the result set iterator */
void kd_res_rewind(struct kdres *set);

/* returns non-zero if the set iterator reached the end after the last element */
int kd_res_end(struct kdres *set);

/* advances the result set iterator, returns non-zero on success, zero if
 * there are no more elements in the result set.
 */
int kd_res_next(struct kdres *set);

/* returns the data pointer (can be null) of the current result set item
 * and optionally sets its position to the pointers(s) if not null.
 */
void *kd_res_item(struct kdres *set, kdcoord *pos);
void *kd_res_item3(struct kdres *set, kdcoord *x, kdcoord *y, kdcoord *z);

/* equivalent to kd_res_item(set, 0) */
void *kd_res_item_data(struct kdres *set);


#ifdef __cplusplus
}
#endif

#endif	/* _KDTREE_H_ */
