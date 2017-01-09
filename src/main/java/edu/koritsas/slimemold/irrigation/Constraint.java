package edu.koritsas.slimemold.irrigation;

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

       Node start =mPath.getFirst();
       Node end =mPath.getLast();

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
                    return dh;
                }
            }));


            SimpleFeature sinkf = (SimpleFeature) end.getObject();
            SimpleFeature sourcef = (SimpleFeature) start.getObject();
            double He= (double) sinkf.getAttribute("hdemand");
        double Ho= (double) sinkf.getAttribute("hdemand");

            double DH= Ho-He;
            //System.out.println("Διαθέσιμο: "+DH+"Πραγματικό: "+dh+" Διαφορά: "+(DH-dh));

             boolean violated=false;
            if ((DH-dh)<0){
                violated=true;

            }
        return violated;
        }

       public boolean containsEdge(Edge edge){

        return mPath.contains(edge);
       }




}
