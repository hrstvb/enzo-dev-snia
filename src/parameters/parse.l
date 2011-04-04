%option noyywrap

D [0-9]
LN [A-Za-z0-9_]
LL [a-z]
LU [A-Z]
LLU [A-Za-z]
E [Ee][+-]?{D}+
S [-]?

%{

/*
 * See LICENSE_CELLO in the main directory for full license agreement
 */

#include "parse.tab.h"

 void yyerror(char *s);

 int    lineno       = 1;


%}

/* double foo (double) */


/* int foo (double) */

/* int    ilogb(double); */
/* int    isnan(double); */

/* double foo (double,double) */

/* atan2     { return ATAN2; } */
/* fmod      { return FMOD; } */
/* hypot     { return HYPOT; } */
/* nextafter { return NEXTAFTER; } */
/* pow       { return POW; } */
/* remainder { return REMAINDER; } */
/* scalb     { return SCALB; } */

/* double jn(int, double); */
/* double ldexp(double, int); */
/* double yn(int, double); */

%%


false     { yylval.logical_type = 0; return LOGICAL; } 
true      { yylval.logical_type = 1; return LOGICAL; } 
[xyzt]    { yylval.string_type = strdup (yytext);
            return VARIABLE; }

{LLU}{LN}* { yylval.string_type = strdup (yytext);
            return IDENTIFIER; }

"="        { return '='; }         /* = */

"*"        { return '*'; }
"/"        { return '/'; }
"+"        { return '+'; }
"-"        { return '-'; }

"<"        { return '<'; }         /* = */
">"        { return '>'; }
"<="        { return LE; }
">="        { return GE; }
"=="        { return EQ; }
"!="        { return NE; }
"&&"        { return AND; }
"||"        { return OR; }
"("          { return '('; }
")"          { return ')'; }

"["        { return '['; } /* begin list */
"]"        { return ']'; } /* end list */
"{"        { return '{'; } /* begin group */
"}"        { return '}'; } /* end group ( and end assignment ) */
";"        { return ';'; } /* parameter assignment separator */
","        { return ','; } /* list element separator */
"#".*      {  /* ignore comments */ }

{S}{D}+{E}            { yylval.scalar_type = atof (yytext); return SCALAR; }
{S}{D}*"."{D}+({E})?  { yylval.scalar_type = atof (yytext); return SCALAR; }
{S}{D}+"."{D}*({E})?  { yylval.scalar_type = atof (yytext); return SCALAR; }
{S}{D}+               { yylval.integer_type = atoi (yytext); return INTEGER; }

\"(\\.|[^\\"])*\"  { yylval.string_type = strdup(yytext+1); /* copy, skipping open quote */
                     yylval.string_type [yyleng-2] = '\0';   /* remove close quote */
                     return STRING; }

[ ] {  }
^[ \t]* {  }

[^ \v\f\t\n]  { /* ignore white space */ }
\n { lineno++; }

. {  }

%%

void yyerror(char *s)
{
  fprintf(stderr, "\n\nERROR: %s in line %d\n\n\n", s,lineno);
  exit(1);
}
