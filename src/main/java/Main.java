import com.vividsolutions.jts.geom.Geometry;

import edu.koritsas.slimemold.TestSlime;
import edu.koritsas.slimemold.mstree.PhysarumPolycephalumLagrarianCSPT;
import edu.koritsas.slimemold.shapefile.GraphUtils;
import edu.koritsas.slimemold.shapefile.IrrigationNetwork;
import edu.koritsas.slimemold.shortestpath.PhysarumPolycephalumLangrarianCSP;
import edu.koritsas.slimemold.shortestpath.PhysarumPolycephalumSP;
import org.geotools.graph.structure.*;
import org.opengis.feature.simple.SimpleFeature;

import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;

/**
 * Created by ilias on 29/9/2016.
 */
public class Main {
    public static void main(String[] args) {

       // IrrigationNetwork network = new IrrigationNetwork("C:/Users/ilias/Desktop/SlimeTest/H.shp","C:/Users/ilias/Desktop/SlimeTest/WS.shp","C:/Users/ilias/Desktop/SlimeTest/P.shp");
        IrrigationNetwork network = new IrrigationNetwork("C:/Users/ilias/Desktop/ParametrizedTests/Hbenchmark1.shp", "C:/Users/ilias/Desktop/ParametrizedTests/Wbenchmark1.shp", "C:/Users/ilias/Desktop/ParametrizedTests/Pbenchmark1.shp");
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
        PhysarumPolycephalumLangrarianCSP slime = new PhysarumPolycephalumLangrarianCSP(graph,source,sinkNodes.get(3),1E-4,1E-4,50000,100,1) {
            @Override
            public boolean pathViolatesConstraints(Graph graph) {
                Node source=getSourceNode();
                Node sink = getSinkNode();

                SimpleFeature sourcef = (SimpleFeature) source.getObject();
                SimpleFeature sinkf = (SimpleFeature) sink.getObject();

                double Ho= (double) sourcef.getAttribute("hdemand");
                double He= (double) sinkf.getAttribute("hdemand");

                double dh = (double) graph.getEdges().stream().collect(Collectors.summingDouble(new ToDoubleFunction<Edge>() {
                    @Override
                    public double applyAsDouble(Edge edge) {


                        return getEdgeConstraintValue(edge);
                    }
                }));
                boolean violates=false;

                if ((Ho-He)<dh){
                    violates=true;
                }

                return violates;
            }

            @Override
            public double getEdgeConstraintValue(Edge edge) {
               SimpleFeature f = (SimpleFeature) edge.getObject();
                double dh = (double) f.getAttribute("Dh");
                return dh;
            }

            @Override
            public double getEdgeCost(Edge e) {
                SimpleFeature f = (SimpleFeature) e.getObject();
                Geometry g = (Geometry) f.getDefaultGeometry();
                double L=g.getLength();
                double cm = (double) f.getAttribute("Cost");
                double cost=L*cm;
                return cost;
            }
        };


        slime.execute();
        slime.showConductivityMap();
        slime.showFlowDiagram();
        GraphUtils.visualizeGraph(graph);
*/

        PhysarumPolycephalumSP slimeSp= new PhysarumPolycephalumSP(graph,source,sinkNodes.get(3),1E-20,1E-20,50000) {
            @Override
            public double getEdgeCost(Edge e) {
                SimpleFeature f = (SimpleFeature) e.getObject();
                Geometry g = (Geometry) f.getDefaultGeometry();
                double L=g.getLength();
                double cm = (double) f.getAttribute("Cost");
                double cost=L*cm;
                return cost;
            }
        };

        slimeSp.execute();
        slimeSp.showFlowDiagram();
        slimeSp.showConductivityMap();
        GraphUtils.visualizeGraph(graph);

    }
}
