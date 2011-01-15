#ifndef MESSAGE_H
#define MESSAGE_H
void c_error (const char *, int );
void c_warning (const char *, int );

#define ERROR_MESSAGE c_error(__FILE__,__LINE__);
#define WARNING_MESSAGE c_warning(__FILE__,__LINE__);
#define DEBUG_MESSAGE c_warning(__FILE__,__LINE__);

      
#endif

