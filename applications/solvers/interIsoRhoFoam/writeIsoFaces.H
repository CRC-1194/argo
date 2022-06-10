void writeIsoFaces
 (
     const fvMesh& mesh_,
     const DynamicList<List<point>>& faces,
     bool writeIsoFacesToFile_
 )
 {
     if (!((writeIsoFacesToFile_ && mesh_.time().writeTime())||(mesh_.time().timeIndex()== 1))) return;
  
     // Writing isofaces to obj file for inspection, e.g. in paraview
     const fileName outputFile
     (
         mesh_.time().globalPath()
       / "isoFaces"
       / word::printf("isoFaces_%012d.obj", mesh_.time().timeIndex())
     );
  
     if (Pstream::parRun())
     {
         // Collect points from all the processors
         List<DynamicList<List<point>>> allProcFaces(Pstream::nProcs());
         allProcFaces[Pstream::myProcNo()] = faces;
         Pstream::gatherList(allProcFaces);
  
         if (Pstream::master())
         {
             mkDir(outputFile.path());
             OBJstream os(outputFile);
             Info<< nl << "isoAdvection: writing iso faces to file: "
                 << os.name() << nl << endl;
  
             face f;
             forAll(allProcFaces, proci)
             {
                 const DynamicList<List<point>>& procFacePts =
                     allProcFaces[proci];
  
                 forAll(procFacePts, i)
                 {
                     const List<point>& facePts = procFacePts[i];
  
                     if (facePts.size() != f.size())
                     {
                         f = face(identity(facePts.size()));
                     }
  
                     os.write(f, facePts, false);
                 }
             }
         }
     }
     else
     {
         mkDir(outputFile.path());
         OBJstream os(outputFile);
         Info<< nl << "isoAdvection: writing iso faces to file: "
             << os.name() << nl << endl;
  
         face f;
         forAll(faces, i)
         {
             const List<point>& facePts = faces[i];
  
             if (facePts.size() != f.size())
             {
                 f = face(identity(facePts.size()));
             }
  
             os.write(f, facePts, false);
         }
     }
 }
