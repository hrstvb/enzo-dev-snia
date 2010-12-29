// $Id: lcamem.h 10 2010-03-18 22:35:19Z jbordner $

/// @file       lcamem.h
/// @author     James Bordner (jobordner@ucsd.edu)
/// @brief      Overrides new() and delete [] () to track memory allocation

#ifndef LCAMEM_H
#define LCAMEM_H

namespace lcamem {
  void *new_(size_t bytes) throw (std::bad_alloc);
  void delete_(void *p) throw ();
  void clear_();
  extern long long bytes_;
  extern long long bytesHigh_;
  extern long long newCalls_;
  extern long long newBytes_;
  extern long long deleteCalls_;
  extern long long deleteBytes_;
}  

extern const bool lcamem_define_arrays_only;

#endif /* LCAMEM_H */
