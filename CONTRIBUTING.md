# 发表指南

## 提交一篇论文

1. 在[总目录](README.md)中选择研究领域。可以复制[论文模板](templates/paper.md)，也可以直接编写普通 Markdown 文件。
2. 将文件放入对应的领域目录，例如 `OPN/my-paper.md`。每篇论文也可以有自己的子目录，和插图、数据或 PDF 附件放在一起。
3. 把改动推送到 `main`，或通过拉取请求合并。拉取请求会先试构建，合并后 GitHub Actions 自动生成领域目录和论文网页。更新、移动或删除稿件也会反映在目录中。

在 [Actions](https://github.com/lixiang90/vibemath-NumberTheory/actions/workflows/pages.yml) 可以查看构建和部署状态。部署成功后，从[在线文库](https://lixiang90.github.io/vibemath-NumberTheory/)进入相应领域即可阅读。

## 标题与署名

文首的 YAML 信息可填写 `title`（标题）、`author`（作者或作者列表）、`date`（发表日期，`YYYY-MM-DD`）和 `abstract`（摘要）。省略这些字段时，系统优先用正文的一级标题作为题名，再以文件名作为后备。日期由作者填写，有日期的稿件按日期从新到旧排列。

每篇新稿默认列入「研究成果」。仅在编写教学或排版样例时添加 `kind: example`，该文将列入「排版示例」。每个领域根目录的 `README.md` 由系统维护，请将论文保存为其他文件名。

## 公式、引文与附件

行内公式写作 `$a^2+b^2=c^2$`，独立公式用两行 `$$` 包住。网站在构建时排好公式，并随站点提供字体。可使用二级标题组织正文，网页侧栏会列出这些章节。

参考文献可列在文末，并附 DOI 或出版社链接。引用仓库内稿件时填写源文件的相对路径，例如 `[相关论文](../OPN/my-paper.md)`。插图写作 `![图注](figures/plot.png)`，附件写作 `[数据](data.csv)`。构建时会将这些链接换成网页或附件地址。

## 添加研究领域

在 [site/topics.json](site/topics.json) 中添加目录名 `id`、中文名 `title` 和简介 `description`，并创建同名目录放入首篇 Markdown。首页、仓库总目录及新领域的成果索引会一起更新。

## 本地预览

安装 Node.js 22 或更高版本，然后在仓库根目录运行以下命令。

```sh
npm ci
npm run build
```

打开生成的 `dist/index.html` 即可阅读。运行 `npm run catalogs` 可在提交前预览 GitHub 上的目录变更。
