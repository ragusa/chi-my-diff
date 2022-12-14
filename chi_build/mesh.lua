--############################################### Setup mesh
chiMeshHandlerCreate()
 
mesh={}
N=4
L=2
xmin = -L/2
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end
 
--chiMeshCreateUnpartitioned3DOrthoMesh(mesh,mesh,mesh)
chiMeshCreateUnpartitioned2DOrthoMesh(mesh,mesh)
--chiMeshCreateUnpartitioned1DOrthoMesh(mesh)
chiVolumeMesherExecute();
 
--############################################### Set Material IDs
--chiVolumeMesherSetMatIDToAll(0)


--############################################### Set Material IDs
material = chiPhysicsAddMaterial("Homogenous_Material");
-- Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,material)

D = {1.1}
Q = {1.2}
XSa = {1.3}
function D_coef(i,x,y,z)
    return D[i+1]
end
function Q_ext(i,x,y,z)
    return Q[i+1]
end
function Sigma_a(i,x,y,z)
    return XSa[i+1]
end


-- Setboundary IDs
-- xmin,xmax,ymin,ymax,zmin,zmax
e_vol = chiLogicalVolumeCreate(RPP,0.99999,1000,-1000,1000,-1000,1000)
w_vol = chiLogicalVolumeCreate(RPP,-1000,-0.9999,-1000,1000,-1000,1000)
n_vol = chiLogicalVolumeCreate(RPP,-1000,1000,0.99999,1000,-1000,1000)
s_vol = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,-0.99999,-1000,1000)

e_bndry = 0
w_bndry = 1
n_bndry = 2
s_bndry = 3

chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,e_vol,e_bndry)
chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,w_vol,w_bndry)
chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,n_vol,n_bndry)
chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,s_vol,s_bndry)

--chiMeshHandlerExportMeshToVTK("Mesh")

--############################################### Add material properties
-- Set material properties
chiPhysicsMaterialAddProperty(material,SCALAR_VALUE,"D")
chiPhysicsMaterialSetProperty(material,"D",SINGLE_VALUE,1.0)

chiPhysicsMaterialAddProperty(material,SCALAR_VALUE,"q")
chiPhysicsMaterialSetProperty(material,"q",SINGLE_VALUE,0.0)

--#### CFEM stuff
phys1 = chiCFEMDiffusionSolverCreate()
chiSolverSetBasicOption(phys1, "residual_tolerance", 1E-8)
chiCFEMDiffusionSetBCProperty(phys1,"boundary_type",e_bndry,"robin", 0.25, 0.5, 0.0)
chiCFEMDiffusionSetBCProperty(phys1,"boundary_type",n_bndry,"reflecting")
chiCFEMDiffusionSetBCProperty(phys1,"boundary_type",s_bndry,"reflecting")
chiCFEMDiffusionSetBCProperty(phys1,"boundary_type",w_bndry,"robin", 0.25, 0.5, 1.0)

-- chiCFEMDiffusionSetProperty(phys1,"boundary_type",n_bndry,"reflecting")
-- chiCFEMDiffusionSetProperty(phys1,"boundary_type",s_bndry,"reflecting")
-- chiCFEMDiffusionSetProperty(phys1,"boundary_type",w_bndry,"neumann",2.0)
-- chiCFEMDiffusionSetProperty(phys1,"boundary_type",e_bndry,"dirichlet",1.0)

chiSolverInitialize(phys1)
chiSolverExecute(phys1)

----############################################### Visualize the field function
fflist,count = chiGetFieldFunctionList(phys1)
chiExportFieldFunctionToVTK(fflist[1],"CFEM","Flux_Diff")