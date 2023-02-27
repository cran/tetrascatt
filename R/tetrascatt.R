  #'tetrascatt
  #'
  #' This function computes the volumetric backscattering from a mesh
  #' of tetrahedrons.

  #' @param parameters a list including  the parameters model, it must include
  #' \itemize{
  #'   \item{\code{cw}: sound speed in the water in m/s}
  #'   \item{\code{g}: g density constrast value, i.e  g= rho1/rhow, where rho1
  #'   and rhow are the density values of the stcatterer and
  #'   the media (sea water) respectively.}
  #'   \item{\code{h}: h density sound speed contrast value, that is
  #'   h= c1/cw where c1 is the sound speed of the stcatterer.}
  #' }
  #'
  #' @param freq  an array of frequencies where the scattering is computed.
  #' @param mesh a list representing the mesh, it must include
  #' \itemize{
  #'   \item{ \code{vertex}: a data.frame with the vertex of the tetrahedra,
  #'   each vertex has to have three coordinates.}
  #'   \item{ \code{tetra}: a data.frame containing the four index
  #'   of each tetrahedra.}
  #' }
  #' @param kversor A three component vector that indicates the direction of the  incident plane wave.
  #'
  #' @return  List containing the frequencies, \code{freq}, and  their corresponding Target Strength values, \code{ts} .
  #' @seealso \code{\link{read_mesh}} to get this kind of list from a .mesh file.
  #' @export
  #' @examples
  #'
  #'
  #' #########################################
  #' ### Set  the Frequency range    ########
  #' #########################################
  #' fmin=12
  #' fmax=400
  #' freqs= seq(fmin,fmax, by=1)
  #' # for tetrascatt freq unities  should be in Hz.
  #' freq=freqs*1000
  #' ############################################################
  #' ########### Set properties of media and scatterer ######
  #' ############################################################
  #' cw <- 1477.4 #soundspeed surrounding fluid (m/s)
  #' rho <- 1026.8 #density surrounding fluid (kg/m^3)
  #' g <- 1028.9/rho #density contrast
  #' h <- 1480.3/cw #soundspeed contrast
  #' my_parameters=list(cw=cw,g=g,h=h)
  #'
  #' ##########################################################
  #' ### Set the incident direction of the plane wave #####
  #' ##########################################################
  #' kversor=c(1,0,0)
  #'
  #' ##########################################################
  #' ### Set the scatterer geometry #######################
  #' ##########################################################
  #' # generates a pseudofile that has  the mesh of cube of one meter
  #' # side
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
  #' # creating an empty temporary mesh file
  #' temp_mesh_file=tempfile(fileext = ".mesh")
  #' # loading the file with data.
  #' writeLines(pseudofile,temp_mesh_file)
  #'
  #'
  #' #reading the mesh
  #' my_mesh=read_mesh( meshfile=temp_mesh_file)
  #'
  #' # Computing the scattering
  #' output= tetrascatt(parameters=my_parameters,freq,
  #'                    mesh=my_mesh,kversor)
  #'
  #' plot(output$freq,output$ts)
  #'
  #' # unliking the teporary file.
  #' unlink(temp_mesh_file)

tetrascatt<-function(parameters,freq,mesh,kversor){
  Tet=mesh$tetra
  Ver=mesh$vertex
  Ver=cbind(Ver[,1],Ver[,2],Ver[,3])
  Tet=cbind(Tet[,1],Tet[,2],Tet[,3],Tet[,4])
  #Tri=my_mesh$triangles
  cw=parameters$cw
  g=parameters$g
  h=parameters$h
#  rhow=parameters$rhow
  finf = tetrascatt_c(cw, g, h, freq, Tet, Ver, kversor )
  ts=20*log10(abs(finf))
  ret=list(freq=freq,ts=ts,finf=finf)
  return(ret)
}

