package edu.koritsas.slimemold.irrigation;

import com.vividsolutions.jts.geom.Geometry;
import org.apache.commons.math3.util.FastMath;
import org.geotools.graph.path.Path;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;
import org.opengis.feature.simple.SimpleFeature;

import java.util.ArrayList;
import java.util.List;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;

/**
 * Created by koritsas on 9/1/2017.
 */
public class Constraint {
    private Path mPath;

    public Constraint(Path branch){
       this.mPath=branch;
    }
    public boolean isViolated(Graph graph){
        List<Edge> edgeList = new ArrayList<>(graph.getEdges());

       Node start =mPath.getLast();
       Node end =mPath.getFirst();

            List<Edge> actual = new ArrayList<>();


            for (Edge e:edgeList){
                if (mPath.contains(e.getNodeA())&&mPath.contains(e.getNodeB())){
                    actual.add(e);
                }
            }


            double dh = actual.stream().collect(Collectors.summingDouble(new ToDoubleFunction<Edge>() {
                @Override
                public double applyAsDouble(Edge edge) {
                       SimpleFeature f = (SimpleFeature) edge.getObject();
                      double dh=(double)f.getAttribute("Dh");
                    double D= (double) f.getAttribute("Diameter");
                    double Q= (double) f.getAttribute("Q");
                    Geometry g = (Geometry) f.getDefaultGeometry();
                    double L=g.getLength();
                    //double dh= 0.00090940294*FastMath.pow(Q,1.78571428571)*L/FastMath.pow(D,4.78571428571);

                    return dh;
                }
            }));


            SimpleFeature sinkf = (SimpleFeature) end.getObject();
            SimpleFeature sourcef = (SimpleFeature) start.getObject();
            double He= (double) sinkf.getAttribute("hdemand");
        double Ho= (double) sourcef.getAttribute("hdemand");

            double DH= Ho-He;
            //System.out.println("Διαθέσιμο: "+DH+"Πραγματικό: "+dh+" Διαφορά: "+(DH-dh));

             boolean violated=false;
            if (DH-dh<-7.8){
                violated=true;

            }
        //System.out.println("Διαθέσιμο: "+DH+"Πραγματικό: "+dh+" Διαφορά: "+(DH-dh));
        return violated;
        }

       public boolean containsEdge(Edge edge){
        boolean contains =false;
        if (mPath.contains(edge.getNodeA())&&mPath.contains(edge.getNodeB())){
            contains=true;
        }


        return contains;
       }
   public double getValue(){
       Node start =mPath.getLast();
       Node end =mPath.getFirst();
       SimpleFeature sinkf = (SimpleFeature) end.getObject();
       SimpleFeature sourcef = (SimpleFeature) start.getObject();
       double He= (double) sinkf.getAttribute("hdemand");
       double Ho= (double) sourcef.getAttribute("hdemand");

       double DH= Ho-He;
       return DH;
   }



}
