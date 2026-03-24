from flask import Flask, render_template, request

app = Flask(__name__)

@app.route("/")
def home():
    return render_template("epigeneticTool.html")

@app.route("/submit", methods=["POST"])
def submit():
    choice = request.form["choice"]
    datafile = request.files["adata"]
    genes = request.form["genes"]
    import anndata as ad
    import scanpy as sc
    import io
    import harmonypy as hm
    import matplotlib.pyplot as plt
    import base64
    from io import BytesIO
    adata = ad.read_h5ad(datafile)
    #make plot based on choice identity
    plot = 8
    adata.obs_names_make_unique()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata) 
    sc.pp.highly_variable_genes(adata, flavor='seurat_v3', layer='X', n_top_genes=3000)
    if choice=="PCA":
        sc.tl.pca(adata,svd_solver='arpack',use_highly_variable=True)
        harmony_result = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'patient_id')
        adata.obsm['X_harmony'] = harmony_result.Z_corr.T 
        sc.pp.neighbors(adata, use_rep='X_harmony', metric='cosine')
        sc.tl.leiden(adata, key_added="leiden_res0_2", resolution=0.2)
        sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
        sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)
        sc.tl.umap(adata)
        sc.pl.umap(adata,color=["leiden_res0_2", "leiden_res0_5", "leiden_res1"],legend_loc="on data",show=False)
        buf = BytesIO()
        plt.savefig(buf, format="png")
        buf.seek(0)
        img = base64.b64encode(buf.read()).decode("utf-8")
        return render_template("result.html", img_data=img)
    elif choice=="t-SNE":
        sc.tl.tsne(adata)
        harmony_result = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'patient_id')
        adata.obsm['X_harmony'] = harmony_result.Z_corr.T 
        sc.pp.neighbors(adata, use_rep='X_harmony', metric='cosine')
        sc.tl.leiden(adata, key_added="leiden_res0_2", resolution=0.2)
        sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
        sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)
        sc.tl.umap(adata)
        sc.pl.umap(adata,color=["leiden_res0_2", "leiden_res0_5", "leiden_res1"],legend_loc="on data",show=False)
        buf = BytesIO()
        plt.savefig(buf, format="png")
        buf.seek(0)
        img = base64.b64encode(buf.read()).decode("utf-8")
        return render_template("result.html", img_data=img)
    elif choice=="heatmap":
        geneList=genes.split(",")
        # sc.tl.pca(adata,svd_solver='arpack',use_highly_variable=True)
        # harmony_result = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'patient_id')
        # adata.obsm['X_harmony'] = harmony_result.Z_corr.T 
        # sc.pp.neighbors(adata, use_rep='X_harmony', metric='cosine')
        # sc.tl.leiden(adata, key_added="leiden_res0_2", resolution=0.2)
        # sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
        # sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)
        # sc.tl.umap(adata)
        sc.pl.heatmap(adata, geneList, groupby='cell_type',show=False)
        buf = BytesIO()
        plt.savefig(buf, format="png")
        buf.seek(0)
        img = base64.b64encode(buf.read()).decode("utf-8")
        return render_template("result.html", img_data=img)
    else:
        sc.pl.correlation_matrix(adata, groupby='cell_type',show=False)
        buf = BytesIO()
        plt.savefig(buf, format="png")
        buf.seek(0)
        img = base64.b64encode(buf.read()).decode("utf-8")
        return render_template("result.html", img_data=img)
if __name__ == "__main__":
    app.run(debug=True)

        #upload epigeneitc data
        #apply visualization algos to it!!
        #What graphing type do you want 
        #UMAP visualization
        #t-SNE visualization
        #Heatmap for marker gene (insert name)
        #HGene-gene correlaiton plots
                # adata
              #what level of processing
              #   $Have python insertion
              #   Have user input group column name, cell type column
              # Do some basic plotting
              #   Recommend plots! Next steps (BIGGEST THINGS)
              #   Cell-level UMAP of adata
          # </py-script>