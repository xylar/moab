#define iMesh_getErrorType FC_FUNC_(imesh_geterrortype, IMESH_GETERRORTYPE)
#define iMesh_getDescription FC_FUNC_(imesh_getdescription, IMESH_GETDESCRIPTION)
#define iMesh_newMesh FC_FUNC_(imesh_newmesh, IMESH_NEWMESH)
#define iMesh_dtor FC_FUNC_(imesh_dtor, IMESH_DTOR)
#define iMesh_load FC_FUNC_(imesh_load, IMESH_LOAD)
#define iMesh_newLoad FC_FUNC_(imesh_newload, IMESH_NEWLOAD)
#define iMesh_save FC_FUNC_(imesh_save, IMESH_SAVE)
#define iMesh_getRootSet FC_FUNC_(imesh_getrootset, IMESH_GETROOTSET)
#define iMesh_getGeometricDimension FC_FUNC_(imesh_getgeometricdimension, IMESH_GETGEOMETRICDIMENSION)
#define iMesh_getDfltStorage FC_FUNC_(imesh_getdfltstorage, IMESH_GETDFLTSTORAGE)
#define iMesh_getAdjTable FC_FUNC_(imesh_getadjtable, IMESH_GETADJTABLE)
#define iMesh_getNumOfType FC_FUNC_(imesh_getnumoftype, IMESH_GETNUMOFTYPE)
#define iMesh_getNumOfTopo FC_FUNC_(imesh_getnumoftopo, IMESH_GETNUMOFTOPO)
#define iMesh_getAllVtxCoords FC_FUNC_(imesh_getallvtxcoords, IMESH_GETALLVTXCOORDS)
#define iMesh_getVtxCoordIndex FC_FUNC_(imesh_getvtxcoordindex, IMESH_GETVTXCOORDINDEX)
#define iMesh_fetErrorType FC_FUNC_(imesh_feterrorType, IMESH_FETERRORTYPE)
#define iMesh_getEntities FC_FUNC_(imesh_getentities, IMESH_GETENTITIES)
#define iMesh_getVtxArrCoords FC_FUNC_(imesh_getvtxarrcoords, IMESH_GETVTXARRCOORDS)
#define iMesh_getAdjEntities FC_FUNC_(imesh_getadjentities, IMESH_GETADJENTITIES)
#define iMesh_initEntArrIter FC_FUNC_(imesh_initentarriter, IMESH_INITENTARRITER)
#define iMesh_getNextEntArrIter FC_FUNC_(imesh_getnextentarriter, IMESH_GETNEXTENTARRITER)
#define iMesh_resetEntArrIter FC_FUNC_(imesh_resetentarriter, IMESH_RESETENTARRITER)
#define iMesh_endEntArrIter FC_FUNC_(imesh_endentarriter, IMESH_ENDENTARRITER)
#define iMesh_getEntArrTopo FC_FUNC_(imesh_getentarrtopo, IMESH_GETENTARRTOPO)
#define iMesh_getEntArrType FC_FUNC_(imesh_getentarrtype, IMESH_GETENTARRTYPE)
#define iMesh_getEntArrAdj FC_FUNC_(imesh_getentarradj, IMESH_GETENTARRADJ)
#define iMesh_createEntSet FC_FUNC_(imesh_createentset, IMESH_CREATEENTSET)
#define iMesh_destroyEntSet FC_FUNC_(imesh_destroyentset, IMESH_DESTROYENTSET)
#define iMesh_isList FC_FUNC_(imesh_islist, IMESH_ISLIST)
#define iMesh_getNumEntSets FC_FUNC_(imesh_getnumentsets, IMESH_GETNUMENTSETS)
#define iMesh_getEntSets FC_FUNC_(imesh_getentsets, IMESH_GETENTSETS)
#define iMesh_addEntToSet FC_FUNC_(imesh_addenttoset, IMESH_ADDENTTOSET)
#define iMesh_rmvEntFromSet FC_FUNC_(imesh_rmventfromset, IMESH_RMVENTFROMSET)
#define iMesh_addEntArrToSet FC_FUNC_(imesh_addentarrtoset, IMESH_ADDENTARRTOSET)
#define iMesh_rmvEntArrFromSet FC_FUNC_(imesh_rmventarrfromset, IMESH_RMVENTARRFROMSET)
#define iMesh_addEntSet FC_FUNC_(imesh_addentset, IMESH_ADDENTSET)
#define iMesh_rmvEntSet FC_FUNC_(imesh_rmventset, IMESH_RMVENTSET)
#define iMesh_isEntContained FC_FUNC_(imesh_isentcontained, IMESH_ISENTCONTAINED)
#define iMesh_isEntSetContained FC_FUNC_(imesh_isentsetcontained, IMESH_ISENTSETCONTAINED)
#define iMesh_addPrntChld FC_FUNC_(imesh_addprntchld, IMESH_ADDPRNTCHLD)
#define iMesh_rmvPrntChld FC_FUNC_(imesh_rmvprntchld, IMESH_RMVPRNTCHLD)
#define iMesh_isChildOf FC_FUNC_(imesh_ischildof, IMESH_ISCHILDOF)
#define iMesh_getNumChld FC_FUNC_(imesh_getnumchld, IMESH_GETNUMCHLD)
#define iMesh_getNumPrnt FC_FUNC_(imesh_getnumprnt, IMESH_GETNUMPRNT)
#define iMesh_getChldn FC_FUNC_(imesh_getchldn, IMESH_GETCHLDN)
#define iMesh_getPrnts FC_FUNC_(imesh_getprnts, IMESH_GETPRNTS)
#define iMesh_setVtxArrCoords FC_FUNC_(imesh_setvtxarrcoords, IMESH_SETVTXARRCOORDS)
#define iMesh_createVtxArr FC_FUNC_(imesh_createvtxarr, IMESH_CREATEVTXARR)
#define iMesh_createEntArr FC_FUNC_(imesh_createentarr, IMESH_CREATEENTARR)
#define iMesh_deleteEntArr FC_FUNC_(imesh_deleteentarr, IMESH_DELETEENTARR)
#define iMesh_createTag FC_FUNC_(imesh_createtag, IMESH_CREATETAG)
#define iMesh_destroyTag FC_FUNC_(imesh_destroytag, IMESH_DESTROYTAG)
#define iMesh_getTagName FC_FUNC_(imesh_gettagname, IMESH_GETTAGNAME)
#define iMesh_getTagSizeValues FC_FUNC_(imesh_gettagsizevalues, IMESH_GETTAGSIZEVALUES)
#define iMesh_getTagSizeBytes FC_FUNC_(imesh_gettagsizebytes, IMESH_GETTAGSIZEBYTES)
#define iMesh_getTagHandle FC_FUNC_(imesh_gettaghandle, IMESH_GETTAGHANDLE)
#define iMesh_getTagType FC_FUNC_(imesh_gettagtype, IMESH_GETTAGTYPE)
#define iMesh_setEntSetData FC_FUNC_(imesh_setentsetdata, IMESH_SETENTSETDATA)
#define iMesh_setEntSetIntData FC_FUNC_(imesh_setentsetintdata, IMESH_SETENTSETINTDATA)
#define iMesh_setEntSetDblData FC_FUNC_(imesh_setentsetdbldata, IMESH_SETENTSETDBLDATA)
#define iMesh_setEntSetBoolData FC_FUNC_(imesh_setentsetbooldata, IMESH_SETENTSETBOOLDATA)
#define iMesh_setEntSetEHData FC_FUNC_(imesh_setentsetehdata, IMESH_SETENTSETEHDATA)
#define iMesh_getEntSetData FC_FUNC_(imesh_getentsetdata, IMESH_GETENTSETDATA)
#define iMesh_getEntSetIntData FC_FUNC_(imesh_getentsetintdata, IMESH_GETENTSETINTDATA)
#define iMesh_getEntSetDblData FC_FUNC_(imesh_getentsetdbldata, IMESH_GETENTSETDBLDATA)
#define iMesh_getEntSetBoolData FC_FUNC_(imesh_getentsetbooldata, IMESH_GETENTSETBOOLDATA)
#define iMesh_getEntSetEHData FC_FUNC_(imesh_getentsetehdata, IMESH_GETENTSETEHDATA)
#define iMesh_getAllEntSetTags FC_FUNC_(imesh_getallentsettags, IMESH_GETALLENTSETTAGS)
#define iMesh_rmvEntSetTag FC_FUNC_(imesh_rmventsettag, IMESH_RMVENTSETTAG)
#define iMesh_setVtxCoords FC_FUNC_(imesh_setvtxcoords, IMESH_SETVTXCOORDS)
#define iMesh_createVtx FC_FUNC_(imesh_createvtx, IMESH_CREATEVTX)
#define iMesh_createEnt FC_FUNC_(imesh_createent, IMESH_CREATEENT)
#define iMesh_deleteEnt FC_FUNC_(imesh_deleteent, IMESH_DELETEENT)
#define iMesh_getArrData FC_FUNC_(imesh_getarrdata, IMESH_GETARRDATA)
#define iMesh_getIntArrData FC_FUNC_(imesh_getintarrdata, IMESH_GETINTARRDATA)
#define iMesh_getDblArrData FC_FUNC_(imesh_getdblarrdata, IMESH_GETDBLARRDATA)
#define iMesh_getBoolArrData FC_FUNC_(imesh_getboolarrdata, IMESH_GETBOOLARRDATA)
#define iMesh_getEHArrData FC_FUNC_(imesh_geteharrdata, IMESH_GETEHARRDATA)
#define iMesh_setArrData FC_FUNC_(imesh_setarrdata, IMESH_SETARRDATA)
#define iMesh_setIntArrData FC_FUNC_(imesh_setintarrdata, IMESH_SETINTARRDATA)
#define iMesh_setDblArrData FC_FUNC_(imesh_setdblarrdata, IMESH_SETDBLARRDATA)
#define iMesh_setBoolArrData FC_FUNC_(imesh_setboolarrdata, IMESH_SETBOOLARRDATA)
#define iMesh_setEHArrData FC_FUNC_(imesh_seteharrdata, IMESH_SETEHARRDATA)
#define iMesh_rmvArrTag FC_FUNC_(imesh_rmvarrtag, IMESH_RMVARRTAG)
#define iMesh_getData FC_FUNC_(imesh_getdata, IMESH_GETDATA)
#define iMesh_getIntData FC_FUNC_(imesh_getintdata, IMESH_GETINTDATA)
#define iMesh_getDblData FC_FUNC_(imesh_getdbldata, IMESH_GETDBLDATA)
#define iMesh_getBoolData FC_FUNC_(imesh_getbooldata, IMESH_GETBOOLDATA)
#define iMesh_getEHData FC_FUNC_(imesh_getehdata, IMESH_GETEHDATA)
#define iMesh_setData FC_FUNC_(imesh_setdata, IMESH_SETDATA)
#define iMesh_setIntData FC_FUNC_(imesh_setintdata, IMESH_SETINTDATA)
#define iMesh_setDblData FC_FUNC_(imesh_setdbldata, IMESH_SETDBLDATA)
#define iMesh_setBoolData FC_FUNC_(imesh_setbooldata, IMESH_SETBOOLDATA)
#define iMesh_setEHData FC_FUNC_(imesh_setehdata, IMESH_SETEHDATA)
#define iMesh_getAllTags FC_FUNC_(imesh_getalltags, IMESH_GETALLTAGS)
#define iMesh_rmvTag FC_FUNC_(imesh_rmvtag, IMESH_RMVTAG)
#define iMesh_initEntIter FC_FUNC_(imesh_initentiter, IMESH_INITENTITER)
#define iMesh_getNextEntIter FC_FUNC_(imesh_getnextentiter, IMESH_GETNEXTENTITER)
#define iMesh_resetEntIter FC_FUNC_(imesh_resetentiter, IMESH_RESETENTITER)
#define iMesh_endEntIter FC_FUNC_(imesh_endentiter, IMESH_ENDENTITER)
#define iMesh_getEntTopo FC_FUNC_(imesh_getenttopo, IMESH_GETENTTOPO)
#define iMesh_getEntType FC_FUNC_(imesh_getenttype, IMESH_GETENTTYPE)
#define iMesh_getVtxCoord FC_FUNC_(imesh_getvtxcoord, IMESH_GETVTXCOORD)
#define iMesh_getEntAdj FC_FUNC_(imesh_getentadj, IMESH_GETENTADJ)
#define iMesh_subtract FC_FUNC_(imesh_subtract, IMESH_SUBTRACT)
#define iMesh_intersect FC_FUNC_(imesh_intersect, IMESH_INTERSECT)
#define iMesh_unite FC_FUNC_(imesh_unite, IMESH_UNITE)
#define iMesh_free FC_FUNC_(imesh_free, IMESH_FREE)
