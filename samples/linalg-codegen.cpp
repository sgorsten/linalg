// The purpose of this test is to check how thoroughly the optimizer is able to inline linalg.h.
// Ideally, no calls should ever exist to linalg:: functions, they should be fully inlined into
// whatever function is making use of them.

// This test is pretty tricky. During compilation, we use volatile variables to force the
// compiler to compile the test functions we want, and to not optimize away the computations
// we wish to evaluate. The program must be compiled with /FAs, to produce a .asm output
// with source annotations. The program then parses its own assembly output, looking for
// specific functions, and counts up the various opcodes found within.

// Finally, we return success or failure on the basis of whether or not each test function
// made the expected number of function calls, which would typically be zero for purely
// algebraic expressions, or some nonzero number for expressions which contain calls to
// sqrt or other math functions.

//////////////
// Preamble //
//////////////

#include <string>
#include <vector>
#include <map>

struct record 
{ 
    static std::vector<record *> records;
    void (*func)(); const char * name; int expected_calls; bool found; std::map<std::string, int> opcodes; std::vector<std::string> calls; 
    record(void (*func)(), const char * name, int expected_calls) : func(func), name(name), expected_calls(expected_calls), found() { records.push_back(this); }
};
std::vector<record *> record::records;

#define TEST(F, C) void F(); record r_##F(&F, #F, C); void F()

///////////
// Tests //
///////////

#include "../linalg.h"
using namespace linalg::aliases;

volatile float v;

TEST(test_dot_float3_float3, 0) // Should fully inline
{
    v = dot(float3(v,v,v), float3(v,v,v)); 
}

TEST(test_nlerp_float3_float3, 1) // Should call sqrtf
{
    auto r = nlerp(float3(v,v,v), float3(v,v,v), v);
    v = r.x; v = r.y; v = r.z;
}

TEST(test_mul_float4x4_float4, 0) // Should fully inline
{
    auto r = linalg::mul(float4x4({v,v,v,v},{v,v,v,v},{v,v,v,v},{v,v,v,v}), float4(v,v,v,v));
    v = r.x; v = r.y; v = r.z; v = r.w;
}

TEST(test_sin_float4, 4) // Should call sinf four times
{
    auto r = sin(float4(v,v,v,v));
    v = r.x; v = r.y; v = r.z; v = r.w;
}

/////////////////
// Test driver //
/////////////////

#include <iostream>
#include <fstream>
#include <sstream>

int main() try
{
    std::string path = "obj/linalg-codegen/Release-";
    path += sizeof(void *) == 8 ? "x64" : "Win32";
    path += "/linalg-codegen.asm";
    std::ifstream in(path);
    if(!in) throw std::runtime_error("Unable to open ./"+path);
    record * current_record = 0;
    std::string line, token0, token1;
    while(true)
    {
        if(!std::getline(in, line)) break;
        if(line.empty()) continue;
        if(line[0] == ';' || line[0] == '$') continue;
        if(line[0] == '?')
        {
            if(!(std::istringstream(line) >> token0 >> token1)) continue;
            if(token1 == "PROC")
            {
                if(current_record) throw std::runtime_error("PROC without matching ENDP");
                auto i = token0.find('@', 1);
                if(i == std::string::npos) throw std::runtime_error("symbol without @");
                auto name = token0.substr(1, i-1);
                for(auto r : record::records)
                {
                    if(name == r->name)
                    {
                        current_record = r;
                        r->found = true;
                        break;
                    }
                }
            }
            else if(token1 == "ENDP")
            {
                current_record = nullptr;
            }
            continue;
        }
        if(!current_record) continue;        
        if(!(std::istringstream(line) >> token0 >> token1)) continue;
        if(token0 == "call" || token0 == "jmp")
        {
            auto i = line.find(';');
            if(i != std::string::npos)
            {
                current_record->calls.push_back(line.substr(i+2));
                continue;
            }
            else if(token0 == "call")
            {
                current_record->calls.push_back(token1);
                continue;
            }
        }
        ++current_record->opcodes[token0];
    }

    int retval = EXIT_SUCCESS;
    void (* volatile v_func)() = 0;
    for(auto * r : record::records)
    {
        if(!r->found) throw std::runtime_error("Function not found: " + std::string(r->name));

        std::cout << r->name << ":\n    ops: ";
        for(auto & p : r->opcodes) std::cout << p.first << "*" << p.second << " ";
        std::cout << "\n    calls: ";
        for(auto & c : r->calls) std::cout << c << " ";
        if(r->calls.empty()) std::cout << "none";
        std::cout << "\n    status: ";
        if(r->calls.size() == r->expected_calls) std::cout << "PASS\n\n";
        else
        {
            retval = EXIT_FAILURE;
            std::cout << "FAIL\n\n";
        }
        v_func = r->func;
    }
    return retval;
}
catch(const std::exception & e)
{
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
}