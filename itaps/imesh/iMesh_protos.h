#ifndef IMESH_PROTOS_H
#define IMESH_PROTOS_H

#include "moab/MOABConfig.h"

#if defined(MOAB_FC_FUNC_)
#define ITAPS_FC_WRAPPER MOAB_FC_FUNC_
#elif defined(MOAB_FC_FUNC)
#define ITAPS_FC_WRAPPER MOAB_FC_FUNC
#else
#define ITAPS_FC_WRAPPER(name,NAME) name
#endif

#define iMesh_getErrorType ITAPS_FC_WRAPPER( imesh_geterrortype, IMESH_GETERRORTYPE )
#define iMesh_getDescription ITAPS_FC_WRAPPER( imesh_getdescription, IMESH_GETDESCRIPTION )
#define iMesh_newMesh ITAPS_FC_WRAPPER( imesh_newmesh, IMESH_NEWMESH )
#define iMesh_dtor ITAPS_FC_WRAPPER( imesh_dtor, IMESH_DTOR )
#define iMesh_load ITAPS_FC_WRAPPER( imesh_load, IMESH_LOAD )
#define iMesh_save ITAPS_FC_WRAPPER( imesh_save, IMESH_SAVE )
#define iMesh_getRootSet ITAPS_FC_WRAPPER( imesh_getrootset, IMESH_GETROOTSET )
#define iMesh_getGeometricDimension ITAPS_FC_WRAPPER( imesh_getgeometricdimension, IMESH_GETGEOMETRICDIMENSION )
#define iMesh_setGeometricDimension ITAPS_FC_WRAPPER( imesh_setgeometricdimension, IMESH_SETGEOMETRICDIMENSION )
#define iMesh_getDfltStorage ITAPS_FC_WRAPPER( imesh_getdfltstorage, IMESH_GETDFLTSTORAGE )
#define iMesh_getAdjTable ITAPS_FC_WRAPPER( imesh_getadjtable, IMESH_GETADJTABLE )
#define iMesh_setAdjTable ITAPS_FC_WRAPPER( imesh_setadjtable, IMESH_SETADJTABLE )
#define iMesh_getNumOfType ITAPS_FC_WRAPPER( imesh_getnumoftype, IMESH_GETNUMOFTYPE )
#define iMesh_getNumOfTopo ITAPS_FC_WRAPPER( imesh_getnumoftopo, IMESH_GETNUMOFTOPO )
#define iMesh_optimize ITAPS_FC_WRAPPER( imesh_optimize, IMESH_OPTIMIZE )
#define iMesh_getEntities ITAPS_FC_WRAPPER( imesh_getentities, IMESH_GETENTITIES )
#define iMesh_getVtxArrCoords ITAPS_FC_WRAPPER( imesh_getvtxarrcoords, IMESH_GETVTXARRCOORDS )
#define iMesh_initEntArrIter ITAPS_FC_WRAPPER( imesh_initentarriter, IMESH_INITENTARRITER )
#define iMesh_getNextEntArrIter ITAPS_FC_WRAPPER( imesh_getnextentarriter, IMESH_GETNEXTENTARRITER )
#define iMesh_resetEntArrIter ITAPS_FC_WRAPPER( imesh_resetentarriter, IMESH_RESETENTARRITER )
#define iMesh_endEntArrIter ITAPS_FC_WRAPPER( imesh_endentarriter, IMESH_ENDENTARRITER )
#define iMesh_getEntArrTopo ITAPS_FC_WRAPPER( imesh_getentarrtopo, IMESH_GETENTARRTOPO )
#define iMesh_getEntArrType ITAPS_FC_WRAPPER( imesh_getentarrtype, IMESH_GETENTARRTYPE )
#define iMesh_getEntArrAdj ITAPS_FC_WRAPPER( imesh_getentarradj, IMESH_GETENTARRADJ )
#define iMesh_getEntArr2ndAdj ITAPS_FC_WRAPPER( imesh_getentarr2ndadj, IMESH_GETENTARR2NDADJ )
#define iMesh_getAdjEntIndices ITAPS_FC_WRAPPER( imesh_getadjentindices, IMESH_GETADJENTINDICES )
#define iMesh_createEntSet ITAPS_FC_WRAPPER( imesh_createentset, IMESH_CREATEENTSET )
#define iMesh_destroyEntSet ITAPS_FC_WRAPPER( imesh_destroyentset, IMESH_DESTROYENTSET )
#define iMesh_isList ITAPS_FC_WRAPPER( imesh_islist, IMESH_ISLIST )
#define iMesh_getNumEntSets ITAPS_FC_WRAPPER( imesh_getnumentsets, IMESH_GETNUMENTSETS )
#define iMesh_getEntSets ITAPS_FC_WRAPPER( imesh_getentsets, IMESH_GETENTSETS )
#define iMesh_addEntToSet ITAPS_FC_WRAPPER( imesh_addenttoset, IMESH_ADDENTTOSET )
#define iMesh_rmvEntFromSet ITAPS_FC_WRAPPER( imesh_rmventfromset, IMESH_RMVENTFROMSET )
#define iMesh_addEntArrToSet ITAPS_FC_WRAPPER( imesh_addentarrtoset, IMESH_ADDENTARRTOSET )
#define iMesh_rmvEntArrFromSet ITAPS_FC_WRAPPER( imesh_rmventarrfromset, IMESH_RMVENTARRFROMSET )
#define iMesh_addEntSet ITAPS_FC_WRAPPER( imesh_addentset, IMESH_ADDENTSET )
#define iMesh_rmvEntSet ITAPS_FC_WRAPPER( imesh_rmventset, IMESH_RMVENTSET )
#define iMesh_isEntContained ITAPS_FC_WRAPPER( imesh_isentcontained, IMESH_ISENTCONTAINED )
#define iMesh_isEntArrContained ITAPS_FC_WRAPPER( imesh_isentarrcontained, IMESH_ISENTARRCONTAINED )
#define iMesh_isEntSetContained ITAPS_FC_WRAPPER( imesh_isentsetcontained, IMESH_ISENTSETCONTAINED )
#define iMesh_addPrntChld ITAPS_FC_WRAPPER( imesh_addprntchld, IMESH_ADDPRNTCHLD )
#define iMesh_rmvPrntChld ITAPS_FC_WRAPPER( imesh_rmvprntchld, IMESH_RMVPRNTCHLD )
#define iMesh_isChildOf ITAPS_FC_WRAPPER( imesh_ischildof, IMESH_ISCHILDOF )
#define iMesh_getNumChld ITAPS_FC_WRAPPER( imesh_getnumchld, IMESH_GETNUMCHLD )
#define iMesh_getNumPrnt ITAPS_FC_WRAPPER( imesh_getnumprnt, IMESH_GETNUMPRNT )
#define iMesh_getChldn ITAPS_FC_WRAPPER( imesh_getchldn, IMESH_GETCHLDN )
#define iMesh_getPrnts ITAPS_FC_WRAPPER( imesh_getprnts, IMESH_GETPRNTS )
#define iMesh_setVtxArrCoords ITAPS_FC_WRAPPER( imesh_setvtxarrcoords, IMESH_SETVTXARRCOORDS )
#define iMesh_createVtxArr ITAPS_FC_WRAPPER( imesh_createvtxarr, IMESH_CREATEVTXARR )
#define iMesh_createEntArr ITAPS_FC_WRAPPER( imesh_createentarr, IMESH_CREATEENTARR )
#define iMesh_deleteEntArr ITAPS_FC_WRAPPER( imesh_deleteentarr, IMESH_DELETEENTARR )
#define iMesh_createTag ITAPS_FC_WRAPPER( imesh_createtag, IMESH_CREATETAG )
#define iMesh_destroyTag ITAPS_FC_WRAPPER( imesh_destroytag, IMESH_DESTROYTAG )
#define iMesh_getTagName ITAPS_FC_WRAPPER( imesh_gettagname, IMESH_GETTAGNAME )
#define iMesh_getTagSizeValues ITAPS_FC_WRAPPER( imesh_gettagsizevalues, IMESH_GETTAGSIZEVALUES )
#define iMesh_getTagSizeBytes ITAPS_FC_WRAPPER( imesh_gettagsizebytes, IMESH_GETTAGSIZEBYTES )
#define iMesh_getTagHandle ITAPS_FC_WRAPPER( imesh_gettaghandle, IMESH_GETTAGHANDLE )
#define iMesh_getTagType ITAPS_FC_WRAPPER( imesh_gettagtype, IMESH_GETTAGTYPE )
#define iMesh_setEntSetData ITAPS_FC_WRAPPER( imesh_setentsetdata, IMESH_SETENTSETDATA )
#define iMesh_setEntSetIntData ITAPS_FC_WRAPPER( imesh_setentsetintdata, IMESH_SETENTSETINTDATA )
#define iMesh_setEntSetDblData ITAPS_FC_WRAPPER( imesh_setentsetdbldata, IMESH_SETENTSETDBLDATA )
#define iMesh_setEntSetEHData ITAPS_FC_WRAPPER( imesh_setentsetehdata, IMESH_SETENTSETEHDATA )
#define iMesh_setEntSetESHData ITAPS_FC_WRAPPER( imesh_setentseteshdata, IMESH_SETENTSETESHDATA )
#define iMesh_getEntSetData ITAPS_FC_WRAPPER( imesh_getentsetdata, IMESH_GETENTSETDATA )
#define iMesh_getEntSetIntData ITAPS_FC_WRAPPER( imesh_getentsetintdata, IMESH_GETENTSETINTDATA )
#define iMesh_getEntSetDblData ITAPS_FC_WRAPPER( imesh_getentsetdbldata, IMESH_GETENTSETDBLDATA )
#define iMesh_getEntSetEHData ITAPS_FC_WRAPPER( imesh_getentsetehdata, IMESH_GETENTSETEHDATA )
#define iMesh_getEntSetESHData ITAPS_FC_WRAPPER( imesh_getentseteshdata, IMESH_GETENTSETESHDATA )
#define iMesh_getAllEntSetTags ITAPS_FC_WRAPPER( imesh_getallentsettags, IMESH_GETALLENTSETTAGS )
#define iMesh_rmvEntSetTag ITAPS_FC_WRAPPER( imesh_rmventsettag, IMESH_RMVENTSETTAG )
#define iMesh_setVtxCoord ITAPS_FC_WRAPPER( imesh_setvtxcoord, IMESH_SETVTXCOORD )
#define iMesh_createVtx ITAPS_FC_WRAPPER( imesh_createvtx, IMESH_CREATEVTX )
#define iMesh_createEnt ITAPS_FC_WRAPPER( imesh_createent, IMESH_CREATEENT )
#define iMesh_deleteEnt ITAPS_FC_WRAPPER( imesh_deleteent, IMESH_DELETEENT )
#define iMesh_getArrData ITAPS_FC_WRAPPER( imesh_getarrdata, IMESH_GETARRDATA )
#define iMesh_getIntArrData ITAPS_FC_WRAPPER( imesh_getintarrdata, IMESH_GETINTARRDATA )
#define iMesh_getDblArrData ITAPS_FC_WRAPPER( imesh_getdblarrdata, IMESH_GETDBLARRDATA )
#define iMesh_getEHArrData ITAPS_FC_WRAPPER( imesh_geteharrdata, IMESH_GETEHARRDATA )
#define iMesh_getESHArrData ITAPS_FC_WRAPPER( imesh_getesharrdata, IMESH_GETESHARRDATA )
#define iMesh_setArrData ITAPS_FC_WRAPPER( imesh_setarrdata, IMESH_SETARRDATA )
#define iMesh_setIntArrData ITAPS_FC_WRAPPER( imesh_setintarrdata, IMESH_SETINTARRDATA )
#define iMesh_setDblArrData ITAPS_FC_WRAPPER( imesh_setdblarrdata, IMESH_SETDBLARRDATA )
#define iMesh_setEHArrData ITAPS_FC_WRAPPER( imesh_seteharrdata, IMESH_SETEHARRDATA )
#define iMesh_setESHArrData ITAPS_FC_WRAPPER( imesh_setesharrdata, IMESH_SETESHARRDATA )
#define iMesh_rmvArrTag ITAPS_FC_WRAPPER( imesh_rmvarrtag, IMESH_RMVARRTAG )
#define iMesh_getData ITAPS_FC_WRAPPER( imesh_getdata, IMESH_GETDATA )
#define iMesh_getIntData ITAPS_FC_WRAPPER( imesh_getintdata, IMESH_GETINTDATA )
#define iMesh_getDblData ITAPS_FC_WRAPPER( imesh_getdbldata, IMESH_GETDBLDATA )
#define iMesh_getEHData ITAPS_FC_WRAPPER( imesh_getehdata, IMESH_GETEHDATA )
#define iMesh_getESHData ITAPS_FC_WRAPPER( imesh_geteshdata, IMESH_GETESHDATA )
#define iMesh_setData ITAPS_FC_WRAPPER( imesh_setdata, IMESH_SETDATA )
#define iMesh_setIntData ITAPS_FC_WRAPPER( imesh_setintdata, IMESH_SETINTDATA )
#define iMesh_setDblData ITAPS_FC_WRAPPER( imesh_setdbldata, IMESH_SETDBLDATA )
#define iMesh_setEHData ITAPS_FC_WRAPPER( imesh_setehdata, IMESH_SETEHDATA )
#define iMesh_setESHData ITAPS_FC_WRAPPER( imesh_seteshdata, IMESH_SETESHDATA )
#define iMesh_getAllTags ITAPS_FC_WRAPPER( imesh_getalltags, IMESH_GETALLTAGS )
#define iMesh_rmvTag ITAPS_FC_WRAPPER( imesh_rmvtag, IMESH_RMVTAG )
#define iMesh_initEntIter ITAPS_FC_WRAPPER( imesh_initentiter, IMESH_INITENTITER )
#define iMesh_getNextEntIter ITAPS_FC_WRAPPER( imesh_getnextentiter, IMESH_GETNEXTENTITER )
#define iMesh_resetEntIter ITAPS_FC_WRAPPER( imesh_resetentiter, IMESH_RESETENTITER )
#define iMesh_endEntIter ITAPS_FC_WRAPPER( imesh_endentiter, IMESH_ENDENTITER )
#define iMesh_getEntTopo ITAPS_FC_WRAPPER( imesh_getenttopo, IMESH_GETENTTOPO )
#define iMesh_getEntType ITAPS_FC_WRAPPER( imesh_getenttype, IMESH_GETENTTYPE )
#define iMesh_getVtxCoord ITAPS_FC_WRAPPER( imesh_getvtxcoord, IMESH_GETVTXCOORD )
#define iMesh_getEntAdj ITAPS_FC_WRAPPER( imesh_getentadj, IMESH_GETENTADJ )
#define iMesh_getEnt2ndAdj ITAPS_FC_WRAPPER( imesh_getent2ndadj, IMESH_GETENT2NDADJ )
#define iMesh_subtract ITAPS_FC_WRAPPER( imesh_subtract, IMESH_SUBTRACT )
#define iMesh_intersect ITAPS_FC_WRAPPER( imesh_intersect, IMESH_INTERSECT )
#define iMesh_unite ITAPS_FC_WRAPPER( imesh_unite, IMESH_UNITE )

#ifndef MOAB_NO_ITAPS_EXTENSIONS

#define iMesh_getEntitiesRec ITAPS_FC_WRAPPER( imesh_getentitiesrec, IMESH_GETENTITIESREC )
#define iMesh_getNumOfTypeRec ITAPS_FC_WRAPPER( imesh_getnumoftyperec, IMESH_GETNUMOFTYPEREC )
#define iMesh_getNumOfTopoRec ITAPS_FC_WRAPPER( imesh_getnumoftoporec, IMESH_GETNUMOFTOPOREC )
#define iMesh_getEntsByTagsRec ITAPS_FC_WRAPPER( imesh_getentsbytagsrec, IMESH_GETENTSBYTAGSREC )
#define iMesh_getEntSetsByTagsRec ITAPS_FC_WRAPPER( imesh_getentsetsbytagsrec, IMESH_GETENTSETSBYTAGSREC )
#define iMesh_MBCNType ITAPS_FC_WRAPPER( imesh_mbcntype, IMESH_MBCNTYPE )
#define iMesh_tagIterate ITAPS_FC_WRAPPER( imesh_tagiterate, IMESH_TAGITERATE )
#define iMesh_connectIterate ITAPS_FC_WRAPPER( imesh_connectiterate, IMESH_CONNECTITERATE )
#define iMesh_coordsIterate ITAPS_FC_WRAPPER( imesh_coordsiterate, IMESH_COORDSITERATE )
#define iMesh_stepEntIter ITAPS_FC_WRAPPER( imesh_stepentiter, IMESH_STEPENTITER )
#define iMesh_stepEntArrIter ITAPS_FC_WRAPPER( imesh_stepentarriter, IMESH_STEPENTARRITER )
#define iMesh_initEntArrIterRec ITAPS_FC_WRAPPER( imesh_initentarriterrec, IMESH_INITENTARRITERREC )
#define iMesh_getAllIfaceTags ITAPS_FC_WRAPPER( imesh_getallifacetags, IMESH_GETALLIFACETAGS )
#define iMesh_createTagWithOptions ITAPS_FC_WRAPPER( imesh_createtagwithoptions, IMESH_CREATETAGWITHOPTIONS )
#define iMesh_createStructuredMesh ITAPS_FC_WRAPPER( imesh_createstructuredmesh, IMESH_CREATESTRUCTUREDMESH )
#define iMesh_freeMemory ITAPS_FC_WRAPPER( imesh_freememory, IMESH_FREEMEMORY )

#endif

#endif
