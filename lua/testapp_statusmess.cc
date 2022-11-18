#include "ChiLua/chi_lua.h"
#include "chi_runtime.h"
#include "chi_log.h"
 
int chiPrintStatus(lua_State *L)
{
    chi::log.Log() << "\nHello from lua function";

    return 0;
}