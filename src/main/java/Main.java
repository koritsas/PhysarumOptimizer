import com.vividsolutions.jts.geom.Geometry;
import edu.koritsas.slimemold.shapefile.GraphUtils;
import edu.koritsas.slimemold.shapefile.IrrigationNetwork;
import org.geotools.graph.build.GraphBuilder;
import org.geotools.graph.build.GraphGenerator;
import org.geotools.graph.build.basic.BasicGraphBuilder;
import org.geotools.graph.build.basic.BasicGraphGenerator;
import org.geotools.graph.build.line.BasicLineGraphGenerator;
import org.geotools.graph.path.DijkstraShortestPathFinder;
import org.geotools.graph.path.Path;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Graphable;
import org.geotools.graph.structure.Node;
import org.geotools.graph.structure.basic.BasicGraph;
import org.geotools.graph.traverse.standard.DijkstraIterator;
import org.opengis.feature.simple.SimpleFeature;

import java.io.IOException;
import java.util.List;

/**
 * Created by ilias on 29/9/2016.
 */
public class Main {
    public static void main(String[] args) {

        IrrigationNetwork network = new IrrigationNetwork("C:/Users/ilias/Desktop/SlimeTest/H.shp","C:/Users/ilias/Desktop/SlimeTest/WS.shp","C:/Users/ilias/Desktop/SlimeTest/P.shp");

        Graph graph=null;
        try {
            graph =network.getGraph();
        } catch (IOException e) {
            e.printStackTrace();
        }

        List<Node> sourceNodes =network.getWaterSource();
        List<Node> sinkNodes = network.getHydrants();
        //sinkNodes=sinkNodes.subList(7,8);

        SlimeMold slimeMold = new SlimeMold(graph,sourceNodes,sinkNodes,2,1.8,50000);
        slimeMold.execute();

        slimeMold.showFlowDiagram();
        slimeMold.showConductivityMap();

        try {
            GraphUtils.visualizeGraph(graph);
        } catch (IOException e) {
            e.printStackTrace();
        }

        System.out.println(slimeMold.getSolutionCost());

/*
        DijkstraIterator.EdgeWeighter weighter = new DijkstraIterator.EdgeWeighter() {
            @Override
            public double getWeight(Edge e) {
                SimpleFeature f = (SimpleFeature) e.getObject();
                Geometry geom = (Geometry) f.getDefaultGeometry();
                return geom.getLength();
            }
        };

        DijkstraShortestPathFinder pf = new DijkstraShortestPathFinder(graph,sourceNodes.get(0),weighter);
        pf.calculate();

        Path path=pf.getPath(sinkNodes.get(0));

        BasicGraph g =new BasicGraph(null,path.getEdges());

        try {
            GraphUtils.visualizeGraph(g);
        } catch (IOException e) {
            e.printStackTrace();
        }

        System.out.println(pf.getCost(sinkNodes.get(0)));
        */
    }
}
