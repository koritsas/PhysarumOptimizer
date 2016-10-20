package edu.koritsas.slimemold.minimumtree;

import org.geotools.graph.build.GraphBuilder;
import org.geotools.graph.build.basic.BasicGraphBuilder;
import org.geotools.graph.path.DijkstraShortestPathFinder;
import org.geotools.graph.path.Path;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;
import org.geotools.graph.structure.basic.BasicGraph;
import org.geotools.graph.traverse.standard.DijkstraIterator;

import java.util.List;

/**
 * Created by ilias on 18/10/2016.
 */
public class DijkstraMinimumTree {
    private Graph graph;
    private Node source;
    private List<Node> destinations;
    private DijkstraIterator iterator;
    private DijkstraIterator.EdgeWeighter edgeWeighter;
    private Graph newGraph;
    private Path path;
    public DijkstraMinimumTree(Graph graph, Node source, List<Node> destinations, DijkstraIterator iterator,DijkstraIterator.EdgeWeighter edgeWeighter){
        this.graph=graph;
        this.source=source;
        this.destinations=destinations;
        this.edgeWeighter=edgeWeighter;
        this.iterator=iterator;
    }
    public void calculate(){

        path=new Path();
        DijkstraShortestPathFinder pf = new DijkstraShortestPathFinder(graph,source,edgeWeighter);
        pf.calculate();
        for (Node d:destinations){
           Path p=pf.getPath(d);
            List<Edge> edges =p.getEdges();
           for (Edge e:edges){
               if (!path.contains(e)){
                   path.getEdges().add(e);

               }
           }

        }
    }
    public Graph getGraph(){
        BasicGraph graph = new BasicGraph(null,path.getEdges());
        return graph;
    }
}
