import com.vividsolutions.jts.geom.Geometry;
import edu.koritsas.slimemold.DirectedSlimeSP;
import edu.koritsas.slimemold.SlimeMoldSP;
import edu.koritsas.slimemold.minimumtree.DijkstraMinimumTree;
import edu.koritsas.slimemold.minimumtree.PhysarumPolycephalumMST;
import edu.koritsas.slimemold.shapefile.DirectedIrrigationNetwork;
import edu.koritsas.slimemold.shapefile.GraphUtils;
import edu.koritsas.slimemold.shapefile.IrrigationNetwork;
import org.geotools.data.collection.ListFeatureCollection;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.graph.path.DijkstraShortestPathFinder;
import org.geotools.graph.path.Path;
import org.geotools.graph.structure.DirectedGraph;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;
import org.geotools.graph.structure.basic.BasicGraph;
import org.geotools.graph.traverse.standard.DijkstraIterator;
import org.geotools.graph.traverse.standard.DirectedDijkstraIterator;
import org.geotools.map.FeatureLayer;
import org.geotools.map.MapContent;
import org.geotools.styling.SLD;
import org.geotools.styling.Style;
import org.geotools.swing.JMapFrame;
import org.opengis.feature.simple.SimpleFeature;

import java.awt.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Created by ilias on 29/9/2016.
 */
public class Main {
    public static void main(String[] args) {

       // IrrigationNetwork network = new IrrigationNetwork("C:/Users/ilias/Desktop/SlimeTest/H.shp","C:/Users/ilias/Desktop/SlimeTest/WS.shp","C:/Users/ilias/Desktop/SlimeTest/P.shp");
       // IrrigationNetwork network = new IrrigationNetwork("C:/Users/ilias/Desktop/ParametrizedTests/H4.shp", "C:/Users/ilias/Desktop/ParametrizedTests/W4.shp", "C:/Users/ilias/Desktop/ParametrizedTests/P4.shp");
        DirectedIrrigationNetwork network = new DirectedIrrigationNetwork("C:/Users/ilias/Desktop/ParametrizedTests/H1.shp", "C:/Users/ilias/Desktop/ParametrizedTests/W1.shp", "C:/Users/ilias/Desktop/ParametrizedTests/P1.shp");
        DirectedGraph graph=null;
        try {
            graph =network.getBasicGraph();
        } catch (IOException e) {
            e.printStackTrace();
        }


        List<Node> sourceNodes =network.getWaterSource();
        List<Node> sinkNodes = network.getHydrants();
        //sinkNodes=sinkNodes.subList(7,8);

        Node source =sourceNodes.get(0);

        Random random = new Random();
        Node sink =sinkNodes.get(random.nextInt(sinkNodes.size()));
/*
        DirectedSlimeSP slimeMold = new DirectedSlimeSP(graph,source,sink,2,1.8,600);
       // SlimeMold slimeMold = new SlimeMold(graph,sourceNodes.get(0),sinkNodes,2,1.8,5000);
        slimeMold.execute();

        slimeMold.showFlowDiagram();
        slimeMold.showConductivityMap();

        try {
            GraphUtils.visualizeGraph(graph);
        } catch (IOException e) {
            e.printStackTrace();
        }

        System.out.println(slimeMold.getSolutionCost());

*/
        DirectedDijkstraIterator.EdgeWeighter weighter = new DirectedDijkstraIterator.EdgeWeighter() {
            @Override
            public double getWeight(Edge e) {
                SimpleFeature f = (SimpleFeature) e.getObject();
                Geometry geom = (Geometry) f.getDefaultGeometry();
                return geom.getLength();
            }
        };
/*
        DijkstraShortestPathFinder pf = new DijkstraShortestPathFinder(graph,source,weighter);
        pf.calculate();

        Path path=pf.getPath(sink);

        BasicGraph g =new BasicGraph(null,path.getEdges());

        try {
            GraphUtils.visualizeGraph(g);
        } catch (IOException e) {
            e.printStackTrace();
        }

        System.out.println(pf.getCost(sink));
*/
        DijkstraMinimumTree tree = new DijkstraMinimumTree(graph,source,sinkNodes,null,weighter);
        tree.calculate();
        Graph dijkstramst =tree.getGraph();



        PhysarumPolycephalumMST slimeTree = new PhysarumPolycephalumMST(graph,source,sinkNodes,2,1.8,500) {
            @Override
            public double getEdgeCost(Edge e) {
                SimpleFeature f = (SimpleFeature) e.getObject();
                Geometry g= (Geometry) f.getDefaultGeometry();
                return g.getLength();
            }
        };

        slimeTree.execute();
        slimeTree.showConductivityMap();
        slimeTree.showFlowDiagram();

        Graph slimeTreeGraph=slimeTree.getGraph();

      GraphUtils.visualizeGraph(slimeTreeGraph);
        GraphUtils.visualizeGraph(dijkstramst);
    }
}
