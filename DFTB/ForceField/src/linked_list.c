/*

 */
#include <stdlib.h>

#include "linked_list.h"

LIST *list_new(void)
{
  LIST *list = (LIST *) malloc(sizeof(LIST));
  list->head = NULL;
  list->tail = NULL;
  return list;
}

void list_append(LIST *list, void *element)
{
  // create a new node
  NODE *node = (NODE *) malloc(sizeof(NODE));
  node->next = NULL;
  node->element = element;
  
  if (list->tail == NULL) {
    list->tail = node;
    list->head = list->tail;
    
  } else {
    list->tail->next = node;
    list->tail = list->tail->next;
  }
}

size_t list_length(LIST *list)
{
  NODE *current;
  size_t len=0;
  for(current=list->head; current!=NULL; current=current->next) {
    len++;
  }
  return len;
}

void list_delete(LIST *list)
{
  NODE *old,*new;
  old = list->head;
  while (old != NULL) {
    new = old->next;
    // free memory of data
    if (old->element) free(old->element);
    // free memory of node
    free(old);
    old = new;
  }
  // free list pointer
  free(list);
}
