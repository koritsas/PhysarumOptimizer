import com.vividsolutions.jts.geom.Geometry;

import edu.koritsas.slimemold.mstree.PhysarumPolycephalumLagrarianCSPT;
import edu.koritsas.slimemold.mstree.PhysarumPolycephalumMST;
import edu.koritsas.slimemold.mstree.PhysarumPolycephalumSPT;
import edu.koritsas.slimemold.shapefile.GraphUtils;
import edu.koritsas.slimemold.shapefile.IrrigationNetwork;
import edu.koritsas.slimemold.shortestpath.PhysarumPolycephalumLangrarianCSP;
import edu.koritsas.slimemold.shortestpath.PhysarumPolycephalumSP;
import org.geotools.graph.structure.*;
import org.opengis.feature.simple.SimpleFeature;

import java.io.IOException;
import java.util.Collection;
import java.util.List;

/**
 * Created by ilias on 29/9/2016.
 */
public class Main {
    public static void main(String[] args) {

       // IrrigationNetwork network = new IrrigationNetwork("C:/Users/ilias/Desktop/SlimeTest/H.shp","C:/Users/ilias/Desktop/SlimeTest/WS.shp","C:/Users/ilias/Desktop/SlimeTest/P.shp");
        IrrigationNetwork network = new IrrigationNetwork("C:/Users/ilias/Desktop/ParametrizedTests/HCSPT.shp", "C:/Users/ilias/Desktop/ParametrizedTests/WCSPT.shp", "C:/Users/ilias/Desktop/ParametrizedTests/PCSPT.shp");
        //DirectedIrrigationNetwork network = new DirectedIrrigationNetwork("C:/Users/ilias/Desktop/ParametrizedTests/HDMST.shp", "C:/Users/ilias/Desktop/ParametrizedTests/WDMST.shp", "C:/Users/ilias/Desktop/ParametrizedTests/PDMST.shp");
       Graph graph=null;
        try {
            graph =network.getBasicGraph();
        } catch (IOException e) {
            e.printStackTrace();
        }


        List<Node> sourceNodes =network.getWaterSource();
        List<Node> sinkNodes = network.getHydrants();
        //sinkNodes=sinkNodes.subList(7,8);

        Node source =sourceNodes.get(0);

/*
        PhysarumPolycephalumSPT slimeSP = new PhysarumPolycephalumSPT(graph,source,sinkNodes,0.00000001,0.000000001,500) {
            @Override
            public double getEdgeCost(Edge e) {
                SimpleFeature f = (SimpleFeature) e.getObject();
                Geometry g = (Geometry) f.getDefaultGeometry();

                return g.getLength();
            }
        };
        slimeSP.execute();
        slimeSP.showConductivityMap();
        slimeSP.showFlowDiagram();
        GraphUtils.visualizeGraph(graph);
*/

/*

 PhysarumPolycephalumLangrarianCSP slimeLag = new PhysarumPolycephalumLangrarianCSP(graph,source,sinkNodes.get(0),0.00000001,0.000000001,50000,10,0.0001) {
     @Override
     public boolean pathViolatesConstraints(Graph graph) {
        double cost=0;
        boolean violates =false;
         Collection<Edge> edges =graph.getEdges();
         for (Edge e:edges){
             cost=cost+getEdgeConstraintValue(e);
         }
         if (cost>100){
             violates=true;
         }

         return violates;
     }

     @Override
     public double getEdgeConstraintValue(Edge edge) {
         SimpleFeature f = (SimpleFeature) edge.getObject();
         double C = (double) f.getAttribute("C");

         return C;
     }

     @Override
     public double getEdgeCost(Edge e) {
         SimpleFeature f = (SimpleFeature) e.getObject();
         Geometry g = (Geometry) f.getDefaultGeometry();

         return g.getLength();
     }
 };

    slimeLag.execute();
    slimeLag.showFlowDiagram();
       slimeLag.showConductivityMap();
        GraphUtils.visualizeGraph(slimeLag.getGraph());
*/


       PhysarumPolycephalumLagrarianCSPT slimeCSPT = new PhysarumPolycephalumLagrarianCSPT(graph,source,sinkNodes,0.00000001,0.000000001,50000,100,0.1) {
            @Override
            public boolean pathViolatesConstraints(Graph graph) {
                double cost=0;
                boolean violates =false;
                Collection<Edge> edges =graph.getEdges();
                for (Edge e:edges){
                    cost=cost+getEdgeConstraintValue(e);
                }
                if (cost>100){
                    violates=true;
                }

                return violates;
            }

            @Override
            public double getEdgeConstraintValue(Edge edge) {
                SimpleFeature f = (SimpleFeature) edge.getObject();
                double C = (double) f.getAttribute("C");

                return C;
            }

            @Override
            public double getEdgeCost(Edge e) {
                SimpleFeature f = (SimpleFeature) e.getObject();
                Geometry g = (Geometry) f.getDefaultGeometry();

                return g.getLength();
            }
        };

        slimeCSPT.execute();
        slimeCSPT.showConductivityMap();
        slimeCSPT.showFlowDiagram();
        GraphUtils.visualizeGraph(slimeCSPT.getGraph());


    }
}
