/*

 */
#include <stdlib.h>

typedef struct NODE
{
  struct NODE *next;
  void *element;
} NODE;

typedef struct LIST
{
  NODE *head;
  NODE *tail;
} LIST;

LIST *list_new(void);
void list_delete(LIST *list);
void list_append(LIST *list, void *element);
size_t list_length(LIST *list);

