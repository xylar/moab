

unix {
  CONFIG += warn_on opengl qt debug
  UI_DIR = .ui
  MOC_DIR = .moc
  OBJECTS_DIR = .obj
}















































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































TEMPLATE	= app
LANGUAGE	= C++

CONFIG	+= qt warn_on release

LIBS	+= -lQVTK -L/home/tjtautg/vtkSNL/bin -lvtksnlIO -lvtksnlGraphics -L/usr/local/ParaViewComplete/VTK/bin -lvtkCommon -lvtkFiltering -lvtkGraphics -lvtkHybrid -lvtkImaging -lvtkIO -lvtkRendering -lvtkfreetype -lvtkftgl -L/home/tjtautg/MOAB -lMOAB -L/cubit/netcdf/lib -lnetcdf_c++ -lnetcdf -L/usr/local/lib/graphviz -ldotneato -lgraph
INCLUDEPATH	+= /home/tjtautg/MOAB /usr/local/ParaViewComplete/VTK/Common /usr/local/ParaViewComplete/VTK/Rendering /usr/local/ParaViewComplete/VTK/Graphics /usr/local/ParaViewComplete/VTK /home/tjtautg/vtkSNL /home/tjtautg/vtkSNL/IO /usr/local/ParaViewComplete/VTK/Filtering /usr/local/include/graphviz

SOURCES	+= main.cpp \
	CropTool.cpp \
	DrawDual.cpp
FORMS	= uiQVDual.ui \
	CropToolpopup.ui

