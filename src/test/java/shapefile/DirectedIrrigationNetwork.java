package shapefile;

import com.vividsolutions.jts.geom.Geometry;
import org.geotools.data.FileDataStore;
import org.geotools.data.FileDataStoreFinder;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.data.simple.SimpleFeatureSource;
import org.geotools.graph.build.feature.FeatureGraphGenerator;
import org.geotools.graph.build.line.BasicDirectedLineGraphBuilder;
import org.geotools.graph.build.line.DirectedLineStringGraphGenerator;
import org.geotools.graph.structure.DirectedGraph;
import org.geotools.graph.structure.GraphVisitor;
import org.geotools.graph.structure.Graphable;
import org.geotools.graph.structure.Node;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Created by ilias on 13/10/2016.
 */
public class DirectedIrrigationNetwork {
    private  String hydrantShapefilePath;
    private  String waterSourceShapefilePath;
    private  String pipeShapefilePath;
    private SimpleFeatureType waterSourceType;
    private SimpleFeatureType hydrantType;
    private SimpleFeatureCollection hydrantFeatureCollection;
    private SimpleFeatureCollection waterSourceFeatureCollection;
    private SimpleFeatureCollection pipeFeatureCollection;

    private DirectedGraph graph;
    public DirectedIrrigationNetwork(String hydrantShapefilePath, String waterSourceShapefilePath, String pipeShapefilePath){
        this.hydrantShapefilePath=hydrantShapefilePath;
        this.pipeShapefilePath=pipeShapefilePath;
        this.waterSourceShapefilePath=waterSourceShapefilePath;

    }
    private void createNetworkFromShp() throws IOException {
        File hydrantFile =new File(hydrantShapefilePath);
        FileDataStore hydrantDataStore = FileDataStoreFinder.getDataStore(hydrantFile);
        SimpleFeatureSource hydrantFeatureSource =hydrantDataStore.getFeatureSource();
        SimpleFeatureCollection hydrants =hydrantFeatureSource.getFeatures();

        hydrantType=hydrantDataStore.getSchema();
        hydrantDataStore.dispose();


        File watersourceFile =new File(waterSourceShapefilePath);
        FileDataStore waterSourceDataStore = FileDataStoreFinder.getDataStore(watersourceFile);
        SimpleFeatureSource waterSourceFeatureSource =waterSourceDataStore.getFeatureSource();
        SimpleFeatureCollection waterSimpleFeatureCollection =waterSourceFeatureSource.getFeatures();
        waterSourceType=waterSourceDataStore.getSchema();
        waterSourceDataStore.dispose();


        File pipeFile = new File(pipeShapefilePath);
        FileDataStore pipeDataStore = FileDataStoreFinder.getDataStore(pipeFile);
        SimpleFeatureSource pipeFeatureSource =pipeDataStore.getFeatureSource();
        SimpleFeatureCollection pipeSimpleFeatureCollection =pipeFeatureSource.getFeatures();
        pipeDataStore.dispose();


        //BasicLineGraphGenerator basicLineGraphGenerator = new BasicLineGraphGenerator();
      /*  LineStringGraphGenerator basicLineGraphGenerator = new LineStringGraphGenerator();
        FeatureGraphGenerator graphGenerator = new FeatureGraphGenerator(basicLineGraphGenerator);
        graphGenerator.setGraphBuilder(new BasicLineGraphBuilder());
*/

        DirectedLineStringGraphGenerator directedLineStringGraphGenerator = new DirectedLineStringGraphGenerator();
        FeatureGraphGenerator graphGenerator = new FeatureGraphGenerator(directedLineStringGraphGenerator);
        graphGenerator.setGraphBuilder(new BasicDirectedLineGraphBuilder());


        SimpleFeatureIterator pipeIt =pipeSimpleFeatureCollection.features();
        while (pipeIt.hasNext()){
            graphGenerator.add(pipeIt.next());

        }
        pipeIt.close();


        SimpleFeatureIterator hit = hydrants.features();
        while(hit.hasNext()){
            SimpleFeature f =hit.next();
            Geometry geom =(Geometry) f.getDefaultGeometry();

            directedLineStringGraphGenerator.getNode(geom.getCoordinate()).setObject(f);


        }
        hit.close();
        SimpleFeatureIterator wsit =waterSimpleFeatureCollection.features();
        while(wsit.hasNext()){
            SimpleFeature f =wsit.next();
            Geometry geom =(Geometry) f.getDefaultGeometry();

            directedLineStringGraphGenerator.getNode(geom.getCoordinate()).setObject(f);


        }
        wsit.close();

        graph = (DirectedGraph) graphGenerator.getGraphBuilder().getGraph();


    }
    public List<Node> getWaterSource(){
        GraphVisitor visitor =new GraphVisitor() {
            public int visit(Graphable component) {
                try {
                    SimpleFeature f = (SimpleFeature) component.getObject();
                    if (f.getType().equals(waterSourceType)) {
                        return 1;
                    } else {
                        return -1;
                    }
                }catch (Exception e){
                    return -1;
                }


            }
        };

        return graph.queryNodes(visitor);
    }

    public List<Node> getHydrants(){
        GraphVisitor visitor =new GraphVisitor() {
            public int visit(Graphable component) {
                Node n = (Node) component;
                try {
                    SimpleFeature f = (SimpleFeature) component.getObject();
                    if (f.getType().equals(hydrantType)) {
                        return 1;
                    } else {
                        return -1;
                    }
                }catch (Exception e){
                    return -1;
                }
            }
        };
        return graph.queryNodes(visitor);
    }



    public DirectedGraph getBasicGraph() throws IOException {
        createNetworkFromShp();

        return graph;
    }


}
