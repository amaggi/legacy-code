#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
        static struct node
	{ int key; struct node *next; };
	static struct node *head, *z, *t;
	void stackinit_() 
	   {
	     head = (struct node *) malloc(sizeof *head);
	     z = (struct node *) malloc(sizeof *z);
	     head->next = z; head->key=0;
	     z->next = z;
	     z->key = 0;
	   }
	void push_(int *p)
	   {
	     int v;
	     v = *p;
	     t = (struct node *) malloc(sizeof *t);	
	     t->key = v; t->next = head->next;	
	     head->next =t;	
	   }
	void pop_(int *x)
	   {
	     t = head->next; head->next = t->next;
	     *x = t->key;
	     free(t);
	   }
	void stackempty_(int *i)
	  { 
	    *i = 0;
            if(head->next == z) *i = 1;
          }
        void stackflush_()
           {
             free(head);
             free(z);
           }
