import fs from "node:fs";
import path from "node:path";
import { collect, encodePath, pages, renderContent, root } from "./content.mjs";
import { articleView, homeView, topicView } from "./views.mjs";

const { topics, guide, assets, routes } = collect();

function markdownLabel(value) {
  return value.replace(/[\\[\]|]/g, "\\$&").replace(/\s*\n\s*/g, " ");
}

function saveCatalog(file, value) {
  const dest = path.join(root, file);
  const previous = fs.existsSync(dest) ? fs.readFileSync(dest, "utf8").replace(/\r\n/g, "\n") : "";
  if (previous !== value) {
    fs.mkdirSync(path.dirname(dest), { recursive: true });
    fs.writeFileSync(dest, value);
  }
}

function links(papers, topic) {
  return papers.map(paper => {
    const source = encodePath(path.posix.relative(topic.id, paper.source));
    const web = `${pages}/${encodePath(paper.output)}`;
    const credit = [paper.author, paper.date].filter(Boolean).map(markdownLabel).join(" · ");
    return `- [${markdownLabel(paper.title)}](<${source}>) · [在线阅读](<${web}>)${credit ? ` · ${credit}` : ""}`;
  }).join("\n");
}

function writeCatalogs() {
  const table = ["| 研究领域 | 在线目录 |", "| --- | --- |"];
  for (const topic of topics) {
    const label = markdownLabel(topic.title);
    const url = `${pages}/${encodePath(topic.id)}/`;
    table.push(`| [${label}](${encodePath(topic.id)}/) | [阅读该领域](${url}) |`);
    const papers = topic.papers.filter(paper => !paper.example);
    const examples = topic.papers.filter(paper => paper.example);
    const sections = [
      `# ${topic.title}`, topic.description,
      `[在线目录](${url}) · [返回总目录](../README.md) · [论文模板](../templates/paper.md)`,
      "## 研究成果", links(papers, topic) || "本领域尚未收录研究成果。"
    ];
    if (examples.length) sections.push("## 排版示例", links(examples, topic));
    saveCatalog(`${topic.id}/README.md`, sections.join("\n\n") + "\n");
  }

  const readme = fs.readFileSync(path.join(root, "README.md"), "utf8").replace(/\r\n/g, "\n");
  const block = /<!-- topics:start -->[\s\S]*?<!-- topics:end -->/;
  if (!block.test(readme)) throw new Error("README.md is missing its topic table markers");
  saveCatalog("README.md", readme.replace(block, () => ["<!-- topics:start -->", ...table, "<!-- topics:end -->"].join("\n")));
}

function buildSite() {
  const dist = path.join(root, "dist");
  fs.rmSync(dist, { recursive: true, force: true });
  function dest(file) {
    const target = path.join(dist, file);
    fs.mkdirSync(path.dirname(target), { recursive: true });
    return target;
  }
  const write = (file, html) => fs.writeFileSync(dest(file), html);
  const copy = (source, file) => fs.copyFileSync(path.join(root, source), dest(file));
  copy("site/site.css", "assets/site.css");
  copy("node_modules/katex/dist/katex.min.css", "assets/katex/katex.min.css");
  fs.cpSync(path.join(root, "node_modules/katex/dist/fonts"), path.join(dist, "assets/katex/fonts"), { recursive: true });
  copy("node_modules/katex/LICENSE", "assets/katex/LICENSE");
  copy("templates/paper.md", "templates/paper.md");
  for (const asset of assets) copy(asset.source, asset.output);
  write("index.html", homeView(topics));
  write(guide.output, articleView(guide, null, renderContent(guide, routes)));
  for (const topic of topics) {
    write(topic.output, topicView(topic));
    for (const paper of topic.papers) {
      write(paper.output, articleView(paper, topic, renderContent(paper, routes)));
    }
  }
}

if (process.argv.includes("--catalogs")) writeCatalogs();
else buildSite();
const count = topics.reduce((n, topic) => n + topic.papers.length, 0);
console.log(`${topics.length} topics, ${count} manuscripts (including examples).`);
