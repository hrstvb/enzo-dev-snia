%{
/*
 * See LICENSE_CELLO in the main directory for full license agreement
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#ifndef __APPLE__
	#include <malloc.h>
#endif
	
#define YYDEBUG 1

  /* Quiet a few -Wall errors */

int yylex (void);
void yyrestart  (FILE * input_file );
void yyerror(char *s);
void yylex_destroy();

#include "parse.tab.h"

#include "parse.h"

const char * node_name[] = {
  "node_unknown",
  "node_operation",
  "node_scalar",
  "node_integer",
  "node_variable",
  "node_function"
  };


const char * op_name[] = {
    "+",
    "-",
    "*",
    "/",
    "<=",
    "<",
    ">=",
    ">",
    "==",
    "!=",
    "&&",
    "||"};

  /* ANY CHANGES HERE MUST BE REFLECTED IN parse.h enum_parameter[] */
  const char * parameter_name[]  = {
    "unknown",
    "sentinel",
    "group",
    "integer",
    "scalar",
    "string",
    "identifier",
    "logical",
    "list",
    "function" };

  /* Structure for storing a single parameter / value pair in a linked list */

  /* The head of the linked list of parameter / value pairs */

  struct param_struct * param_head = NULL; /* head of entire list */
  struct param_struct * param_curr = NULL; /* head of current list */

  /* The current groups and parameter type */


  char *              current_parameter = NULL;
  char *              current_group[MAX_GROUP_DEPTH];
  int                 current_group_level = 0;
  enum enum_parameter current_type      = enum_parameter_sentinel;

  void clear_groups (char * groups[]) {
    int i;
    for (i=0; i<MAX_GROUP_DEPTH; i++) {
      groups[i] = 0; 
    }
  };

  void copy_groups (char * group_dest[], char * group_src[]) {
    int i;
    for (i=0; i<MAX_GROUP_DEPTH; i++) {
      group_dest[i] = (group_src[i]) ? strdup(group_src[i]) : 0;
    }
  };

  /* Function to update parameter's groups once the group is known */

/*   void update_group (char * group) */
/*     { */
/*       struct param_struct * p = param_curr; */
/*       while (p->next->type  != enum_parameter_sentinel &&  */
/* 	     p->next->group == NULL) { */
/* 	p->next->group = strdup(group); */
/*         p = p -> next; */
/*       } */
/*     } */

  /* Insert a parameter into the list */

  void insert_param(struct param_struct * head, struct param_struct * new)
  {
     new->next  = head->next;
     head->next = new;
  }

  /* Delete a parameter from the list given a pointer to the previous element */

  void delete_param(struct param_struct * previous)
  {
    struct param_struct * item = previous->next;
    previous->next = item->next;
    free (item);     
  }

  /* Function to update parameter's subgroups once the subgroup is known */

/*   void update_subgroup (char * subgroup) */
/*     { */
/*       struct param_struct * p = param_curr; */
/*       int inside_subgroup = 1; */
/*       while (p->next->type     != enum_parameter_sentinel &&  */
/* 	     p->next->subgroup == NULL) { */
/* 	if (p->next->type == enum_parameter_subgroup) { */
/* 	  inside_subgroup = 0; */
/*           delete_param(p); */
/*         } else if (inside_subgroup) { */
/*           p->next->subgroup = strdup(subgroup); */
/*           p = p -> next; */
/*         } */
/*       } */
/*     } */

  struct param_struct * reverse_param(struct param_struct * old_head)
  {
    /* Keep sentinel the same */

    struct param_struct * new_head = old_head;

    struct param_struct * p = old_head;
    struct param_struct * c = p->next;
    struct param_struct * n = c->next;

    do {
      /* If parameter is a list, recursively reverse it as well */
      if (c->type == enum_parameter_list) {
	c->list_value = reverse_param(c->list_value);
      }
      c->next = p;
      p       = c;
      c       = n;
      n = n->next;
    } while (p->type != enum_parameter_sentinel) ;

    new_head = p;
    return new_head;
  }

  /* Function for creating and inserting a new parameter / value pair */
  /* in the linked list */

  struct param_struct * new_param ()
  {
    /* Create the new node */

    /* MEMORY LEAK */
     struct param_struct * p = 
       (struct param_struct *) malloc (sizeof (struct param_struct));

   /* Fill in the non-type-specific values for the new node */

     /* MEMORY LEAK */
     
     copy_groups(p->group,current_group);

     p->parameter = (current_parameter) ? strdup(current_parameter) : 0;

     current_type = enum_parameter_unknown;

     insert_param(param_curr,p);

   /* Clear variables for the next assignment */

     return p;
  }

  /* New integer parameter assignment */

  void new_param_integer (int value)
  {
    struct param_struct * p = new_param();
    p->type          = enum_parameter_integer;
    p->integer_value = value;
  }


  /* New scalar parameter assignment */

  void new_param_scalar (double value)
  {
    struct param_struct * p = new_param();
    p->type         = enum_parameter_scalar;
    p->scalar_value = value;
  }

  /* New logical parameter assignment */

  void new_param_logical (int value)
  {
    struct param_struct * p = new_param();
    p->type          = enum_parameter_logical;
    p->logical_value = value;
  }

  /* New string parameter assignment */
  void new_param_string (char * value)
  {
    struct param_struct * p = new_param();
    p->type         = enum_parameter_string;
    p->string_value = value;
  }

  /* New subgroup  */
  void new_param_group (char * value)
  {
    struct param_struct * p = new_param();
    p->type         = enum_parameter_group;
    p->string_value = value;
  }

  /* New empty parameter assignment: FIRST NODE IN LIST IS A SENTINEL  */
  struct param_struct * new_param_sentinel ()
  {

    /* MEMORY LEAK */
    struct param_struct * p = 
      (struct param_struct *) malloc (sizeof (struct param_struct));

    clear_groups(p->group);
    p->parameter = NULL;
    p->type      = enum_parameter_sentinel;
    p->next       = p;
    p->list_value = NULL;

    return p;
  }

  /* New list parameter assignment */

  void new_param_list (struct param_struct * curr)
  {
    struct param_struct * p = new_param();
    p->type       = enum_parameter_list;
    p->list_value = curr;
  }

  void new_parameter()
  {
     switch (current_type) {
     case enum_parameter_group:
       new_param_group(yylval.group_type);
       break;
     case enum_parameter_integer:
       new_param_integer(yylval.integer_type);
       break;
     case enum_parameter_scalar:
       new_param_scalar(yylval.scalar_type);
       break;
     case enum_parameter_string: 
       new_param_string(yylval.string_type);
       break;
     case enum_parameter_logical:
       new_param_logical(yylval.logical_type);
       break;
     case enum_parameter_list:
       break;
    default:
       printf ("%s:%d Parse Error: unknown type %d\n",
	       __FILE__,__LINE__,current_type);
       break;
     }
  }

  char * strcat3 (const char * s1,const char * s2,const char * s3)
  {
    char * s = malloc (strlen(s1) + strlen(s2) + strlen(s3) + 1);

    strcpy(s,s1);
    strcpy(s+strlen(s1),s2);
    strcpy(s+strlen(s1)+strlen(s2),s3);
    return s;
  }

  char * ftoa (double f)
    { 
      char * a = malloc(25); 

      sprintf (a,"%24.16e",f);
      return a;
    }

%}


%union { 
  int logical_type;  
  int integer_type; 
  double scalar_type;  
  char * string_type; 
  char * group_type;
  }

/* %token <string_type>  GROUP_NAME */
%token <string_type>  STRING
%token <string_type>  IDENTIFIER
%token <string_type> VARIABLE
%token <scalar_type>  SCALAR
%token <integer_type> INTEGER
%token <logical_type> LOGICAL

%token LE
%token GE
%token NE
%token EQ
%token AND
%token OR

%left OR
%left AND
%left EQ NE 
%left LE GE '<' '>'
%left '+' '-'
%left '*' '/'

/* double foo (double) */

/* double foo (double,double) */

 /* %token ATAN2 FMOD HYPOT NEXTAFTER POW REMAINDER SCALB */

%%

file : /* nothing */
 | file group { }
 ;

group: 
group_name parameter_group  {  }

parameter_group :
'{' parameter_list '}'          { current_group[--current_group_level] = 0; }
| '{' parameter_list ';' '}'    { current_group[--current_group_level] = 0; }

parameter_list : 
                      parameter_assignment  {  }
 | parameter_list ';' parameter_assignment  {  }
 |                    group  {  }
 | parameter_list ';' group  {  }

 
group_name :
  IDENTIFIER                    { current_group[current_group_level++] = $1; }

parameter_name :
  IDENTIFIER                    { current_parameter = $1;} 

parameter_assignment : 
  parameter_name '=' parameter_value { new_parameter(); }
 ;

parameter_value : 
 STRING { current_type = enum_parameter_string;       yylval.string_type = $1; }
 | INTEGER  { current_type = enum_parameter_integer;      yylval.integer_type = $1;}
 | SCALAR   { current_type = enum_parameter_scalar;       yylval.scalar_type = $1;}
 | LOGICAL  { current_type = enum_parameter_logical;      yylval.logical_type = $1; }
 | list { current_type = enum_parameter_list; }
 ;

list: LIST_BEGIN list_elements LIST_END {  }
    | LIST_BEGIN LIST_END {  }

LIST_BEGIN:
 '[' { 
   struct param_struct * p = new_param_sentinel();
   p->list_value = param_curr;
   new_param_list(p);
   param_curr = p;
 }
LIST_END:
 ']' { param_curr = param_curr->list_value; }


list_elements:
                      parameter_value    { new_parameter(); }
 | list_elements  ',' parameter_value    { new_parameter(); }

{ }
;

%%

struct param_struct * 
cello_parameters_read(FILE * fp)
{
  clear_groups(current_group);

  /* initialize the linked list with an initial sentinel (sentinel) node */
  param_head = param_curr = new_param_sentinel();

  /*   yydebug=1; */
  
  yyrestart(fp);

  yyparse();
  yylex_destroy();

  param_head = reverse_param(param_head);
  return param_head;
}

void indent (int level)
{
  int i;
  for (i=0; i<level; i++) {
    printf ("  "); 
  }
}

void cello_parameters_print_list(struct param_struct * head, int level)
{
  struct param_struct * p = head->next;
  int i;

  while (p && p->type != enum_parameter_sentinel) {

    if (p->group != NULL) {
      indent(level);
      printf ("%s ", parameter_name[p->type]);
      for (i=0; p->group[i] != NULL && i < MAX_GROUP_DEPTH; i++) {
	printf ("%s:",p->group[i]);
      }
      printf ("%s = ", p->parameter);
    } else {
      /* list element */
      indent(level);
      printf ("%s %s = ", 
	      parameter_name[p->type], p->parameter);
    }
    switch (p->type) {
    case enum_parameter_scalar:  
      printf ("%g\n",p->scalar_value);  
      break;
    case enum_parameter_integer: 
      printf ("%d\n",p->integer_value); 
      break;
    case enum_parameter_string:  
      printf ("%s\n",p->string_value); 
      break;
    case enum_parameter_group:  
      printf ("Uh oh: GROUP %s (should be deleted)\n",p->string_value);
      break;
    case enum_parameter_logical:
      printf ("%s\n",p->logical_value ? "true" : "false");
      break;
    case enum_parameter_list:    
      indent(level);
      printf ("[\n"); 
      cello_parameters_print_list(p->list_value, level + 1);
      indent(level);
      printf ("]\n"); 
      break;
    default: 
      indent(level);
      printf ("unknown type\n"); 
      break;
    }
    p = p->next;
  }
}

void cello_parameters_print()
{
  cello_parameters_print_list(param_head,0);
}

