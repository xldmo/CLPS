/* Deterministic public test inputs; not keys for operational use. */
#include <stdio.h>
#include <string.h>
#include "clps.h"
#include "params.h"

#define REQUIRE_OK(call) do { if ((call)!=0) { \
    fprintf(stderr,"Failed at line %d: %s\n",__LINE__,#call); return 1; } } while(0)

int main(void) {
    uint8_t master[SEEDBYTES], seedp[SEEDBYTES], seeda[SEEDBYTES];
    uint8_t params[CLPS_PARAMSBYTES], msk[CLPS_MSKBYTES];
    uint8_t partp[CLPS_PARTIALKEYBYTES], parta[CLPS_PARTIALKEYBYTES];
    uint8_t pkp[CLPS_PKIDBYTES], pka[CLPS_PKIDBYTES];
    uint8_t skp[CLPS_SKIDBYTES], ska[CLPS_SKIDBYTES];
    uint8_t token[CLPS_PROXYKEYBYTES], sig[CLPS_SIGBYTES];
    uint8_t idp[SEEDBYTES], ida[SEEDBYTES], uskp[SEEDBYTES], uska[SEEDBYTES];
    uint8_t message[32];
    for(unsigned i=0;i<SEEDBYTES;i++) { master[i]=(uint8_t)i; seedp[i]=(uint8_t)(i^0x11); seeda[i]=(uint8_t)(i^0x22); }
    for(unsigned i=0;i<sizeof(message);i++) message[i]=(uint8_t)i;
    memset(idp,1,sizeof(idp)); memset(ida,2,sizeof(ida));
    memset(uskp,3,sizeof(uskp)); memset(uska,4,sizeof(uska));
    REQUIRE_OK(clps_setup(params,msk,master));
    REQUIRE_OK(clps_partial_key_gen(partp,params,msk,idp));
    REQUIRE_OK(clps_partial_key_gen(parta,params,msk,ida));
    REQUIRE_OK(clps_full_key_gen(pkp,skp,params,partp,idp,uskp,seedp));
    REQUIRE_OK(clps_full_key_gen(pka,ska,params,parta,ida,uska,seeda));
    REQUIRE_OK(clps_principal_authorize_agent(token,params,skp,idp,ida,pka));
    REQUIRE_OK(clps_agent_sign(sig,message,sizeof(message),params,token,pkp,ska));
    REQUIRE_OK(clps_agent_verify(message,sizeof(message),sig,params,pkp,pka,ida,idp));
    message[0]^=1;
    if(clps_agent_verify(message,sizeof(message),sig,params,pkp,pka,ida,idp)==0) { fputs("Altered message accepted\n",stderr); return 1; }
    message[0]^=1; ida[0]^=1;
    if(clps_agent_verify(message,sizeof(message),sig,params,pkp,pka,ida,idp)==0) { fputs("Altered identity accepted\n",stderr); return 1; }
    ida[0]^=1; sig[0]^=1;
    if(clps_agent_verify(message,sizeof(message),sig,params,pkp,pka,ida,idp)==0) { fputs("Altered signature accepted\n",stderr); return 1; }
    printf("CLPS mode %d: valid workflow and three rejection checks passed\n",DILITHIUM_MODE);
    return 0;
}
