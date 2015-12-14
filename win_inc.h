// Include windows.h without evil and overload stuff

#pragma once

#define WIN32_LEAN_AND_MEAN
#define NOATOM
#define NOGDI
#define NOGDICAPMASKS
#define NOMETAFILE
#define NOMINMAX
#define NOMSG
#define NOOPENFILE
#define NORASTEROPS
#define NOSCROLL
#define NOSOUND
#define NOSYSMETRICS
#define NOTEXTMETRIC
#define NOWH
#define NOCOMM
#define NOKANJI
#define NOCRYPT
#define NOMCX
//#define UNICODE
#define STRICT

#include <windows.h>

#undef WIN32_LEAN_AND_MEAN
#undef NOATOM
#undef NOGDI
#undef NOGDICAPMASKS
#undef NOMETAFILE
#undef NOMINMAX
#undef NOMSG
#undef NOOPENFILE
#undef NORASTEROPS
#undef NOSCROLL
#undef NOSOUND
#undef NOSYSMETRICS
#undef NOTEXTMETRIC
#undef NOWH
#undef NOCOMM
#undef NOKANJI
#undef NOCRYPT
#undef NOMCX
//#undef UNICODE
#undef STRICT