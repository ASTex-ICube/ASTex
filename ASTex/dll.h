#ifndef __ASTEX_DLL_H_
#define __ASTEX_DLL_H_

/**
* \brief Linkage declaration for ASTex symbols.
*/
#ifdef WIN32 
#ifdef BUILD_STATIC
#define ASTEX_API
#endif
#ifndef ASTEX_API
#if defined ASTEX_DLL_EXPORT
#define ASTEX_API __declspec(dllexport)
#else
#define ASTEX_API __declspec(dllimport)
#endif
#endif
#else
#define ASTEX_API
#endif

#endif // __ASTEX_DLL_H_
