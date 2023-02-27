  #'read_mesh
  #'
  #' This function reads the mesh from a file .mesh (extension)
  #'

  #' @param meshfile a string with the name of the file that contains the
  #' volumetric mesh in GMF format  (Gamma Mesh Format), conventionally,
  #' an ASCII file with ".mesh" extension.
  #'
  #' @return a list representing the mesh, it should include
  #' \itemize{
  #'   \item{ \code{vertex}: a data frame with the vertices of the tetrahedra,
  #'   each vertex must have three coordinates}
  #'   \item{ \code{tetra}: a data frame containing the four vertex-index of
  #'   each tetrahedron}
  #' }
  #' @export
  #' @examples
  #'
  #'
  #'
  #' # Generates a pseudofile that has  the mesh of
  #' # a cube with edges one metre in length, centered at the origin.
  #'
  #' pseudofile=c("MeshVersionFormatted 2",
  #'   "","Dimension 3","","Vertices","8","-0.5 -0.5 0.5 6 ",
  #'    "-0.5 -0.5 -0.5 7 ","-0.5 0.5 0.5 9 ","-0.5 0.5 -0.5 11 ",
  #'    "0.5 -0.5 0.5 16 ","0.5 -0.5 -0.5 17 ","0.5 0.5 0.5 19 ",
  #'    "0.5 0.5 -0.5 21 ","","Edges", "12","2 1 5 ","1 3 8 ",
  #'    "4 3 10 ","2 4 12 ","6 5 15 ","5 7 18 ","8 7 20 ","6 8 22 ",
  #'   "2 6 25 ","1 5 26 ","4 8 29 ","3 7 30 ","","Triangles",
  #'   "12","2 1 3 3 ","3 4 2 3 ","5 6 8 13 ","8 7 5 13 ",
  #'   "2 6 5 23 ","5 1 2 23 ","8 4 3 27 ","3 7 8 27 ","2 4 8 31 ",
  #'   "8 6 2 31 ","3 1 5 33 ","5 7 3 33 ","","Tetrahedra",
  #'   "5","5 2 1 3 1 ","4 2 8 3 1 ","8 5 7 3 1 ","8 2 6 5 1 ",
  #'   "3 2 8 5 1 ","","End","" )
  #'
  #'# creating an empty temporary mesh file
  #'temp_mesh_file=tempfile(fileext = ".mesh")
  #'# loading the file with data.
  #'writeLines(pseudofile,temp_mesh_file)
  #'# reading the mesh
  #'my_mesh=read_mesh( meshfile=temp_mesh_file)
  #'
  #'# see the bounding box of the volumetric mesh.
  #'lapply(my_mesh$vertex,range)
  #'
  #' # unliking the teporary file.
  #' unlink(temp_mesh_file)
  #'
  read_mesh=function(meshfile){
    mesh=utils::read.csv2(meshfile)
    skip=3
    #leo el numero de vertices
    nvert=as.numeric(mesh[skip,1])
    #separo las filas de los vertices
    nv0=skip+1
    vertices=mesh[seq(nv0,nv0+nvert-1),]
    vertices <- strsplit(vertices, " +")
    vertices <- data.frame(t(sapply(vertices,c)))

    #leo el numero de edges
    nedges=as.numeric(mesh[nv0+nvert+1,1])
    #separo los edges
    ne0=nv0+nvert+2
    edges=mesh[seq(ne0,ne0+nedges-1),]
    edges <- strsplit(edges, " +")
    edges <- data.frame(t(sapply(edges,c)))

    #leo el numero de triangles
    ntriangle=as.numeric(mesh[ne0+nedges+1,1])
    #separo los triangles
    nt0=ne0+nedges+2
    triangles=mesh[seq(nt0,nt0+ntriangle-1),]
    triangles <- strsplit(triangles, " +")
    triangles <- data.frame(t(sapply(triangles,c)))

    #leo el numero de tetrahedros
    ntetra=as.numeric(mesh[nt0+ntriangle+1,1])
    #separo los tetrahedros
    nte0=nt0+ntriangle+2
    tetra=mesh[seq(nte0,nte0+ntetra-1),]
    tetra <- strsplit(tetra, " +")
    tetra <- data.frame(t(sapply(tetra,c)))
    vertices <- lapply(vertices,as.numeric)
    tetra <- lapply(tetra,as.numeric)
    names(tetra)<-c("vertex1","vertex2","vertex3","vertex4", "idk")
    #tetra=tetra[,1:4]
    names(vertices) <- c("x","y","z","idk")
    #vertices = vertices[,1:3];
    mesh=list(vertex =as.data.frame(vertices),tetra = as.data.frame(tetra),
              triangles=as.data.frame(triangles))
    mesh
  }
